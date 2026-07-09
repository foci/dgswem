#!/usr/bin/env python3
"""TUI wrapper for adcprep — DG-SWEM domain decomposition pre-processor.

Run from the directory containing fort.14, fort.15, and fort.dg.

Usage:
    python adcprep_tui.py [--adcprep /path/to/build/adcprep]

If --adcprep is omitted the script searches ../build/adcprep,
../build-parallel/adcprep, ./adcprep, and PATH.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import threading
from dataclasses import dataclass, field
from pathlib import Path

from textual import on
from textual.app import App, ComposeResult
from textual.containers import Horizontal, ScrollableContainer, Vertical
from textual.screen import Screen
from textual.widgets import (
    Button,
    Footer,
    Header,
    Input,
    Label,
    Log,
    RadioButton,
    RadioSet,
    Rule,
    Static,
)


# ── Locate adcprep binary ────────────────────────────────────────────────────

def find_adcprep(hint: str | None) -> Path:
    if hint:
        p = Path(hint)
        if p.is_file() and p.stat().st_mode & 0o111:
            return p.resolve()
        sys.exit(f"adcprep not found or not executable at: {hint}")
    for candidate in [
        Path("adcprep"),
        Path("../build/adcprep"),
        Path("../build-parallel/adcprep"),
    ]:
        if candidate.is_file() and candidate.stat().st_mode & 0o111:
            return candidate.resolve()
    which = shutil.which("adcprep")
    if which:
        return Path(which)
    sys.exit(
        "adcprep not found. Build it first, then pass the path:\n"
        "  python adcprep_tui.py --adcprep ../build/adcprep"
    )


# ── Fort.15 / Fort.14 parsers ────────────────────────────────────────────────

def _ival(lines: list[str], idx: int, pos: int = 0) -> int:
    if idx >= len(lines):
        return 0
    parts = lines[idx].split()
    try:
        return int(parts[pos]) if pos < len(parts) else 0
    except (ValueError, IndexError):
        return 0


def parse_fort15(path: Path) -> dict:
    """Extract flags from fort.15 that control which optional files are needed.

    Fort.15 is positional: each line is parsed left-to-right, trailing text
    ignored. All line indices below are 0-based.
    """
    result = dict(nolifa=0, nwp=0, nws=0, nrs=0, ntif=0, nbfr=0, ok=False)
    try:
        lines = path.read_text(errors="replace").splitlines()

        result["nolifa"] = _ival(lines, 9)   # line 10: NOLIFA
        nwp = _ival(lines, 12)               # line 13: NWP
        result["nwp"] = nwp

        # line 16+NWP: raw combined NWS/NRS token
        raw = _ival(lines, 15 + nwp)
        hundreds = abs(raw) // 100
        nrs = hundreds if hundreds in (1, 2, 3) else 0
        nws = (abs(raw) - hundreds * 100) * (1 if raw >= 0 else -1)
        result["nws"] = nws
        result["nrs"] = nrs

        # W: number of extra lines inserted after NWS for wind/wave time specs
        abs_nws = abs(nws)
        if abs_nws in {2, 4, 5, 6, 8, 9, 12, 15, 16, 19, 20, 45}:
            w = 1
        elif abs_nws == 3:
            w = 2
        elif abs_nws == 0 and nrs >= 1:
            w = 1
        else:
            w = 0

        ntif = _ival(lines, 30 + nwp + w)               # line 31+NWP+W: NTIF
        result["ntif"] = ntif
        result["nbfr"] = _ival(lines, 31 + nwp + w + 2 * ntif)  # NBFR
        result["ok"] = True
    except Exception:
        pass
    return result


def parse_fort14_nope(path: Path) -> int:
    """Return NOPE (open boundary segment count) from fort.14."""
    try:
        lines = path.read_text(errors="replace").splitlines()
        ne = int(lines[1].split()[0])
        np_ = int(lines[1].split()[1])
        return int(lines[2 + np_ + ne].split()[0])
    except Exception:
        return 0


# ── Configuration state ──────────────────────────────────────────────────────

@dataclass
class PrepConfig:
    mnproc: int = 1
    job: str = "prep"
    ihots: int = 67
    grid_file: str = "fort.14"
    runinfo_file: str = "fort.15"
    fort12_file: str = "fort.12"
    fort13_file: str = "fort.13"
    fort19_file: str = "fort.19"
    fort20_file: str = "skip"
    fort22_file: str = "fort.22"
    fort23_file: str = "fort.23"
    ft15: dict = field(default_factory=dict)
    nope: int = 0


def build_stdin_lines(cfg: PrepConfig) -> list[str]:
    """Return the exact lines to pipe to adcprep stdin, in order."""
    lines = [str(cfg.mnproc), cfg.job]
    if cfg.job == "hotprep1":
        lines.append(str(cfg.ihots))
        return lines
    lines += [cfg.grid_file, cfg.runinfo_file]
    ft15 = cfg.ft15
    if ft15.get("nolifa") == 3:
        lines.append(cfg.fort12_file)
    lines.append(cfg.fort13_file)
    if cfg.nope > 0 and ft15.get("nbfr", 0) == 0:
        lines.append(cfg.fort19_file)
    if ft15.get("nws", 0) != 0:
        lines.append(cfg.fort22_file)
    if ft15.get("nrs", 0) == 1:
        lines.append(cfg.fort23_file)
    return lines


# ── Mascot ───────────────────────────────────────────────────────────────────

MASCOT = """
   ____
  /.   \\__        Namo here —
 /_   \\_/ \\       can I help you run a coastal flooding simulation?
// \\ ____ |\\      v0.1
  __|_| |_|snd"""


class MascotBanner(Static):
    DEFAULT_CSS = """
    MascotBanner {
        width: 100%;
        height: auto;
        padding: 1 3;
        color: $primary;
        background: $panel;
        border-bottom: wide $surface-lighten-2;
    }
    """

    def __init__(self) -> None:
        super().__init__(MASCOT)


# ── Wizard base screen ────────────────────────────────────────────────────────

class WizardScreen(Screen):
    """Base class that renders the mascot header on every wizard step."""

    CSS = ""  # subclasses set their own CSS

    def compose(self) -> ComposeResult:
        yield Header()
        yield MascotBanner()
        with Vertical(classes="screen-body"):
            yield from self._body()
        yield Footer()

    def _body(self) -> ComposeResult:
        yield from ()


# ── Shared CSS ───────────────────────────────────────────────────────────────

COMMON_CSS = """
Screen {
    layout: vertical;
}

.screen-body {
    width: 100%;
    height: 1fr;
    align: center middle;
}

.wizard-box {
    width: 78;
    max-height: 40;
    border: round $primary;
    padding: 1 2;
    background: $surface;
}

.step-label {
    color: $text-muted;
    text-style: italic;
}

.step-title {
    text-style: bold;
    padding-bottom: 1;
}

.description {
    color: $text-muted;
    padding-bottom: 1;
}

.field-label {
    text-style: bold;
    margin-top: 1;
}

.field-desc {
    color: $text-muted;
}

.nav {
    margin-top: 2;
    align-horizontal: right;
    height: auto;
}

.nav Button {
    margin-left: 1;
}

.error-msg {
    color: $error;
    padding-top: 1;
}

.summary-row {
    padding: 0 0 0 2;
}

.summary-key {
    color: $text-muted;
    width: 28;
}

.summary-val {
    color: $success;
    text-style: bold;
}
"""


# ── Helper ───────────────────────────────────────────────────────────────────

def _file_ok(path_str: str, skippable: bool = False) -> bool:
    if skippable and path_str.strip().lower() == "skip":
        return True
    return Path(path_str).is_file()


# ── Wizard screens ───────────────────────────────────────────────────────────

class ProcessorsScreen(WizardScreen):
    """Step 1: number of processors."""

    CSS = COMMON_CSS

    def _body(self) -> ComposeResult:
        with ScrollableContainer(classes="wizard-box"):
            yield Static("Step 1 of 6 — Processor count", classes="step-label")
            yield Static("Number of parallel processors", classes="step-title")
            yield Static(
                "adcprep partitions the mesh into N subdomains — one per MPI rank.\n"
                "This must match the -np value you will use with mpirun.",
                classes="description",
            )
            yield Label("Number of processors:", classes="field-label")
            yield Input(
                value="4",
                placeholder="e.g. 32",
                id="mnproc",
            )
            yield Static("", id="err", classes="error-msg")
            with Horizontal(classes="nav"):
                yield Button("Next →", id="next", variant="primary")

    @on(Button.Pressed, "#next")
    def handle_next(self) -> None:
        val = self.query_one("#mnproc", Input).value.strip()
        try:
            n = int(val)
            if n < 1:
                raise ValueError
        except ValueError:
            self.query_one("#err", Static).update("Enter a positive integer.")
            return
        self.app.config.mnproc = n
        self.app.push_screen(JobTypeScreen())


class JobTypeScreen(WizardScreen):
    """Step 2: prep / hotprep1 / hotprep2."""

    CSS = COMMON_CSS

    JOB_DESCRIPTIONS = {
        "prep": (
            "Full prep (cold start)",
            "Partition the mesh and distribute all input files to PE* subdirectories.\n"
            "Use this for a new run starting from rest.",
        ),
        "hotprep1": (
            "Hot start — update fort.15 only (hotprep1)",
            "Assumes PE* subdirectories already exist from a previous run segment.\n"
            "Rewrites each local fort.15 to set NWS negative for continued wind forcing.\n"
            "Fastest option when restarting with the same mesh.",
        ),
        "hotprep2": (
            "Hot start — full re-decomposition (hotprep2)",
            "Runs a full prep AND distributes the hot-start file (fort.67 or fort.68)\n"
            "to each PE* subdirectory. Use when resuming after changing processor count.",
        ),
    }

    def _body(self) -> ComposeResult:
        with ScrollableContainer(classes="wizard-box"):
            yield Static("Step 2 of 6 — Job type", classes="step-label")
            yield Static("What kind of decomposition?", classes="step-title")
            with RadioSet(id="job-radio"):
                for key, (label, _) in self.JOB_DESCRIPTIONS.items():
                    yield RadioButton(label, id=f"job-{key}", value=(key == "prep"))
            yield Static("", id="job-desc", classes="description")
            with Horizontal(classes="nav"):
                yield Button("← Back", id="back")
                yield Button("Next →", id="next", variant="primary")

    def on_mount(self) -> None:
        self._update_desc("prep")

    def _update_desc(self, key: str) -> None:
        _, desc = self.JOB_DESCRIPTIONS[key]
        self.query_one("#job-desc", Static).update(desc)

    @on(RadioSet.Changed, "#job-radio")
    def handle_radio(self, event: RadioSet.Changed) -> None:
        key = event.pressed.id.removeprefix("job-")
        self._update_desc(key)

    @on(Button.Pressed, "#back")
    def go_back(self) -> None:
        self.app.pop_screen()

    @on(Button.Pressed, "#next")
    def handle_next(self) -> None:
        radio = self.query_one("#job-radio", RadioSet)
        key = radio.pressed_button.id.removeprefix("job-")
        self.app.config.job = key
        if key == "hotprep1":
            self.app.push_screen(HotRestartScreen())
        else:
            self.app.push_screen(GridFileScreen())


class HotRestartScreen(WizardScreen):
    """Step 3 (hotprep1 branch): choose fort.67 or fort.68."""

    CSS = COMMON_CSS

    def _body(self) -> ComposeResult:
        with ScrollableContainer(classes="wizard-box"):
            yield Static("Step 3 of 3 — Hot-start file unit", classes="step-label")
            yield Static("Which restart file to use?", classes="step-title")
            yield Static(
                "ADCIRC writes hot-start snapshots alternately to fort.67 and fort.68.\n"
                "Choose whichever was written last (check file modification times).",
                classes="description",
            )
            with RadioSet(id="ihots-radio"):
                yield RadioButton("fort.67 (unit 67)", id="ihots-67", value=True)
                yield RadioButton("fort.68 (unit 68)", id="ihots-68")
            with Horizontal(classes="nav"):
                yield Button("← Back", id="back")
                yield Button("Run adcprep →", id="next", variant="primary")

    @on(Button.Pressed, "#back")
    def go_back(self) -> None:
        self.app.pop_screen()

    @on(Button.Pressed, "#next")
    def handle_next(self) -> None:
        radio = self.query_one("#ihots-radio", RadioSet)
        self.app.config.ihots = int(radio.pressed_button.id.removeprefix("ihots-"))
        self.app.push_screen(ConfirmScreen())


class GridFileScreen(WizardScreen):
    """Step 3: fort.14 mesh file path."""

    CSS = COMMON_CSS

    def _body(self) -> ComposeResult:
        with ScrollableContainer(classes="wizard-box"):
            yield Static("Step 3 of 6 — Mesh file", classes="step-label")
            yield Static("ADCIRC Unit 14 (grid file)", classes="step-title")
            yield Static(
                "The finite-element mesh: node coordinates, element connectivity,\n"
                "and boundary segment definitions. Usually named fort.14.",
                classes="description",
            )
            yield Label("Path to fort.14:", classes="field-label")
            yield Input(value="fort.14", id="grid-input")
            yield Static("", id="err", classes="error-msg")
            with Horizontal(classes="nav"):
                yield Button("← Back", id="back")
                yield Button("Next →", id="next", variant="primary")

    @on(Button.Pressed, "#back")
    def go_back(self) -> None:
        self.app.pop_screen()

    @on(Button.Pressed, "#next")
    def handle_next(self) -> None:
        val = self.query_one("#grid-input", Input).value.strip()
        if not Path(val).is_file():
            self.query_one("#err", Static).update(f"File not found: {val}")
            return
        self.app.config.grid_file = val
        self.app.push_screen(RunInfoScreen())


class RunInfoScreen(WizardScreen):
    """Step 4: fort.15 run-info file path. Triggers fort.15 + fort.14 parse."""

    CSS = COMMON_CSS

    def _body(self) -> ComposeResult:
        with ScrollableContainer(classes="wizard-box"):
            yield Static("Step 4 of 6 — Run parameters file", classes="step-label")
            yield Static("ADCIRC Unit 15 (run info file)", classes="step-title")
            yield Static(
                "Controls all run parameters: time step, tidal forcing, wind model,\n"
                "wetting/drying, output requests, and more. Usually named fort.15.",
                classes="description",
            )
            yield Label("Path to fort.15:", classes="field-label")
            yield Input(value="fort.15", id="runinfo-input")
            yield Static("", id="err", classes="error-msg")
            with Horizontal(classes="nav"):
                yield Button("← Back", id="back")
                yield Button("Next →", id="next", variant="primary")

    @on(Button.Pressed, "#back")
    def go_back(self) -> None:
        self.app.pop_screen()

    @on(Button.Pressed, "#next")
    def handle_next(self) -> None:
        val = self.query_one("#runinfo-input", Input).value.strip()
        if not Path(val).is_file():
            self.query_one("#err", Static).update(f"File not found: {val}")
            return
        self.app.config.runinfo_file = val
        # Parse both files to determine which optional prompts are needed
        self.app.config.ft15 = parse_fort15(Path(val))
        self.app.config.nope = parse_fort14_nope(Path(self.app.config.grid_file))
        self.app.push_screen(OptionalFilesScreen())


class OptionalFilesScreen(WizardScreen):
    """Step 5: optional/conditional input files."""

    CSS = COMMON_CSS

    def _build_specs(self) -> list[dict]:
        cfg = self.app.config
        ft15 = cfg.ft15
        specs = []

        # fort.13 — always prompted, always skippable
        specs.append(dict(
            id="fort13",
            label="Nodal Attributes file (fort.13)",
            desc="Spatially-varying Manning's n, surface canopy, primitive weighting.\n"
                 "Many idealized test cases do not use one.",
            default=cfg.fort13_file,
            skippable=True,
        ))

        # fort.12 — only if NOLIFA=3
        if ft15.get("nolifa") == 3:
            specs.append(dict(
                id="fort12",
                label="StartDry file (fort.12)",
                desc="Required when NOLIFA=3: lists nodes that begin in the dry state.",
                default=cfg.fort12_file,
                skippable=False,
            ))

        # fort.19 — if open boundaries exist and no tidal BCs defined
        nope = cfg.nope
        nbfr = ft15.get("nbfr", 0)
        if nope > 0 and nbfr == 0:
            specs.append(dict(
                id="fort19",
                label="Aperiodic Elevation BC file (fort.19)",
                desc=f"Time-series water-level forcing at open boundaries (NOPE={nope}, no tidal constituents).",
                default=cfg.fort19_file,
                skippable=False,
            ))

        # fort.22 — if wind forcing is active
        nws = ft15.get("nws", 0)
        if nws != 0:
            specs.append(dict(
                id="fort22",
                label="Wind / Pressure Stress file (fort.22)",
                desc=f"Meteorological forcing file (NWS={nws}). Format varies by NWS value; see fort.15 docs.",
                default=cfg.fort22_file,
                skippable=False,
            ))

        # fort.23 — if wave radiation stress is active
        if ft15.get("nrs", 0) == 1:
            specs.append(dict(
                id="fort23",
                label="Wave Radiation Stress file (fort.23)",
                desc="Wave radiation stress gradients from a coupled wave model (NRS=1).",
                default=cfg.fort23_file,
                skippable=False,
            ))

        return specs

    CSS = COMMON_CSS + """
    .input-row {
        height: auto;
    }
    .input-row Input {
        width: 1fr;
    }
    .input-row Button {
        width: auto;
        min-width: 10;
        margin-left: 1;
    }
    """

    def _body(self) -> ComposeResult:
        specs = self._build_specs()
        self._specs = specs
        # Preserve pre-skip value so Include restores to what the user typed.
        self._saved: dict[str, str] = {s["id"]: s["default"] for s in specs}

        with ScrollableContainer(classes="wizard-box"):
            yield Static("Step 5 of 6 — Optional input files", classes="step-label")
            yield Static("Ancillary files detected from fort.15", classes="step-title")
            if not specs:
                yield Static("No optional files required for this configuration.", classes="description")
            else:
                yield Static("Enter file paths below.", classes="description")
                for spec in specs:
                    yield Label(spec["label"], classes="field-label")
                    yield Static(spec["desc"], classes="field-desc")
                    if spec["skippable"]:
                        with Horizontal(classes="input-row"):
                            yield Input(value=spec["default"], id=f"opt-{spec['id']}")
                            yield Button("Skip", id=f"skip-{spec['id']}", variant="warning")
                    else:
                        yield Input(value=spec["default"], id=f"opt-{spec['id']}")
            yield Static("", id="err", classes="error-msg")
            yield Rule()
            with Horizontal(classes="nav"):
                yield Button("← Back", id="back")
                yield Button("Review →", id="next", variant="primary")

    @on(Button.Pressed)
    def handle_button(self, event: Button.Pressed) -> None:
        btn_id = event.button.id or ""
        if btn_id == "back":
            self.app.pop_screen()
        elif btn_id == "next":
            self._do_next()
        elif btn_id.startswith("skip-"):
            self._toggle_skip(btn_id.removeprefix("skip-"))

    def _toggle_skip(self, field_id: str) -> None:
        inp = self.query_one(f"#opt-{field_id}", Input)
        btn = self.query_one(f"#skip-{field_id}", Button)
        if inp.disabled:
            inp.disabled = False
            inp.value = self._saved[field_id]
            btn.label = "Skip"
            btn.variant = "warning"
        else:
            self._saved[field_id] = inp.value
            inp.disabled = True
            btn.label = "Include"
            btn.variant = "default"

    def _do_next(self) -> None:
        cfg = self.app.config
        err = self.query_one("#err", Static)

        for spec in self._specs:
            inp = self.query_one(f"#opt-{spec['id']}", Input)
            if inp.disabled:
                val = "skip"
            else:
                val = inp.value.strip()
                if not Path(val).is_file():
                    err.update(f"File not found: {val}")
                    return
            setattr(cfg, f"{spec['id']}_file", val)

        err.update("")
        self.app.push_screen(ConfirmScreen())


class ConfirmScreen(WizardScreen):
    """Step 6: review all choices before running."""

    CSS = COMMON_CSS

    def _body(self) -> ComposeResult:
        cfg = self.app.config
        ft15 = cfg.ft15

        rows: list[tuple[str, str]] = [("Processors", str(cfg.mnproc)), ("Job type", cfg.job)]

        if cfg.job == "hotprep1":
            rows.append(("Hot-start unit", str(cfg.ihots)))
        else:
            rows += [
                ("Mesh (fort.14)", cfg.grid_file),
                ("Run info (fort.15)", cfg.runinfo_file),
            ]
            if ft15.get("ok"):
                rows += [
                    ("  NOLIFA", str(ft15.get("nolifa", "?"))),
                    ("  NWS", str(ft15.get("nws", "?"))),
                    ("  NRS", str(ft15.get("nrs", "?"))),
                    ("  NBFR", str(ft15.get("nbfr", "?"))),
                    ("  NOPE (from fort.14)", str(cfg.nope)),
                ]
            else:
                rows.append(("  fort.15 parse", "failed — defaults used"))

            rows.append(("Nodal attrs (fort.13)", cfg.fort13_file))
            if ft15.get("nolifa") == 3:
                rows.append(("StartDry (fort.12)", cfg.fort12_file))
            if cfg.nope > 0 and ft15.get("nbfr", 0) == 0:
                rows.append(("Elev BCs (fort.19)", cfg.fort19_file))
            if ft15.get("nws", 0) != 0:
                rows.append(("Wind (fort.22)", cfg.fort22_file))
            if ft15.get("nrs", 0) == 1:
                rows.append(("Wave stress (fort.23)", cfg.fort23_file))

        with ScrollableContainer(classes="wizard-box"):
            yield Static("Step 6 of 6 — Review & run", classes="step-label")
            yield Static("Confirm these settings before launching adcprep", classes="step-title")
            yield Rule()
            for key, val in rows:
                with Horizontal(classes="summary-row"):
                    yield Static(key, classes="summary-key")
                    yield Static(val, classes="summary-val")
            yield Rule()
            yield Static(
                f"adcprep: {self.app.adcprep_path}",
                classes="description",
            )
            with Horizontal(classes="nav"):
                yield Button("← Back", id="back")
                yield Button("Run adcprep", id="run", variant="success")

    @on(Button.Pressed, "#back")
    def go_back(self) -> None:
        self.app.pop_screen()

    @on(Button.Pressed, "#run")
    def handle_run(self) -> None:
        self.app.push_screen(RunScreen())


class RunScreen(WizardScreen):
    """Final screen: streams adcprep output live."""

    CSS = COMMON_CSS + """
    RunScreen .wizard-box {
        width: 90;
        max-height: 38;
    }
    Log {
        height: 24;
        border: round $surface-lighten-2;
    }
    """

    def _body(self) -> ComposeResult:
        with Vertical(classes="wizard-box"):
            yield Static("Running adcprep", classes="step-title")
            yield Static("", id="status", classes="description")
            yield Log(id="run-log", highlight=True, auto_scroll=True)
            with Horizontal(classes="nav"):
                yield Button("Close", id="close", variant="primary", disabled=True)

    def on_mount(self) -> None:
        self._start_adcprep()

    def _start_adcprep(self) -> None:
        cfg = self.app.config
        stdin_lines = build_stdin_lines(cfg)
        stdin_text = "\n".join(stdin_lines) + "\n"

        log = self.query_one("#run-log", Log)
        log.write_line(f"$ {self.app.adcprep_path}")
        for answer in stdin_lines:
            log.write_line(f"  > {answer}")
        log.write_line("")

        def run_in_thread():
            try:
                proc = subprocess.Popen(
                    [str(self.app.adcprep_path)],
                    stdin=subprocess.PIPE,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                )
            except Exception as exc:
                self.app.call_from_thread(self._on_done, f"Failed to start adcprep: {exc}", False)
                return

            # Write stdin in a background thread so a fast exit or early
            # EOF from adcprep does not block or kill the stdout reader.
            def _write_stdin():
                try:
                    proc.stdin.write(stdin_text)
                    proc.stdin.flush()
                except (BrokenPipeError, OSError):
                    pass
                finally:
                    try:
                        proc.stdin.close()
                    except OSError:
                        pass

            threading.Thread(target=_write_stdin, daemon=True).start()

            try:
                for line in proc.stdout:
                    self.app.call_from_thread(log.write_line, line.rstrip())
            except Exception:
                pass

            returncode = proc.wait()
            if returncode == 0:
                msg = "adcprep completed successfully."
            else:
                msg = f"adcprep exited with code {returncode}."
            self.app.call_from_thread(self._on_done, msg, returncode == 0)

        threading.Thread(target=run_in_thread, daemon=True).start()

    def _on_done(self, message: str, success: bool) -> None:
        self.query_one("#status", Static).update(message)
        btn = self.query_one("#close", Button)
        btn.disabled = False
        btn.variant = "success" if success else "error"

    @on(Button.Pressed, "#close")
    def handle_close(self) -> None:
        self.app.exit()


# ── Application ──────────────────────────────────────────────────────────────

class AdcprepTUI(App):
    TITLE = "adcprep — DG-SWEM Domain Decomposition"
    SUB_TITLE = "Parallel pre-processor"
    BINDINGS = [("q", "quit", "Quit"), ("escape", "pop_screen", "Back")]

    def __init__(self, adcprep_path: Path):
        super().__init__()
        self.adcprep_path = adcprep_path
        self.config = PrepConfig()

    def on_mount(self) -> None:
        self.push_screen(ProcessorsScreen())


# ── Entry point ──────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Interactive TUI for adcprep domain decomposition."
    )
    parser.add_argument(
        "--adcprep",
        metavar="PATH",
        help="Path to the adcprep executable (default: auto-detect)",
    )
    args = parser.parse_args()
    adcprep_path = find_adcprep(args.adcprep)
    AdcprepTUI(adcprep_path).run()


if __name__ == "__main__":
    main()
