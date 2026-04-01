#!/usr/bin/env python3
import pytest
from util import run_serial, run_parallel
import os


def test_quarter_annular(binpath):
    run_serial(binpath, "quarter_annular", 0.05, 1e-7)

def test_quarter_annular_parallel(binpath):
    run_parallel(binpath, "quarter_annular", 0.05, 1e-7)

def test_wetdry(binpath):
    run_serial(binpath, "wetdry", 0.01, 1e-7)

def test_wetdry_parallel(binpath):
    run_parallel(binpath, "wetdry", 0.01, 1e-7)

def test_mass_conservation(binpath):
    run_serial(binpath, "mass_conservation", 1e-3, 1e-7)
