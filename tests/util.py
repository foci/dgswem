#!/usr/bin/env python3

import numpy as np
import pandas as pd
import sys
import subprocess
import pytest
import os


def mean_rel_err(xs, ys):
    ab = np.abs(xs - ys)
    mask = ab < 90000
    e = ab[mask]
    e = e / np.maximum(np.abs(xs[mask]), np.abs(ys[mask]))
    return np.mean(e)

def max_abs_err(xs, ys):
    e = np.abs(xs - ys)
    e = e[e < 90000]
    return np.max(e)


def mean_abs_err(xs, ys):
    e = np.abs(xs - ys)
    # Remove wet/dry outliers
    e = e[e < 90000]
    return np.mean(e)

def last_snapshot(filename):
    f = pd.read_csv(filename,
                    sep=r"\s+",
                    lineterminator='\n',
                    skiprows=0)
    num_nodes = int(f.iloc[0,1])
    last = f.iloc[-num_nodes:, 1].to_numpy()
    return last, num_nodes

def run_serial(binpath, testpath, rtol=0.05, atol=0.01):
    if not os.path.exists(os.path.join(binpath, "dgswem_serial")):
        raise FileNotFoundError('dgswem_serial not found')

    result = subprocess.run(os.path.join(binpath, "dgswem_serial"),
                            check=False, cwd=testpath, shell=True,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT)
    if result.returncode != 0:
        print(result.stdout.decode())
        sys.exit("DG-SWEM serial run failed. See output above.")

    d1, _ = last_snapshot(os.path.join(testpath , "fort.63.true"))
    d2, _ = last_snapshot(os.path.join(testpath , "fort.63"))
    np.testing.assert_allclose(d2, d1, rtol=rtol, atol=atol)

def run_parallel(binpath, testpath, rtol=0.05, atol=0.01):
    if not os.path.exists(os.path.join(binpath, "dgswem")):
        pytest.skip("dgswem executable not found. Skipping...")
    if not os.path.exists(os.path.join(binpath, "adcprep")):
        pytest.skip("adcprep executable not found. Skipping...")
    if not os.path.exists(os.path.join(binpath, "adcpost")):
        pytest.skip("adcpost executable not found. Skipping...")

    result = subprocess.run(os.path.join(binpath, "adcprep") + " < in.prep", check=False, cwd=testpath, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        print(result.stdout.decode())
        sys.exit("ADCPREP failed. See output above.")


    result = subprocess.run("mpirun -np 2 " + os.path.join(binpath, "dgswem"), check=False, cwd=testpath, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        print(result.stdout.decode())
        sys.exit("DG-SWEM parallel run failed. See output above.")

    result = subprocess.run(os.path.join(binpath, "adcpost") + " < out.prep", check=False, cwd=testpath, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        print(result.stdout.decode())
        sys.exit("ADCPOST failed. See output above.")

    d1, _ = last_snapshot(os.path.join(testpath , "fort.63.true"))
    d2, _ = last_snapshot(os.path.join(testpath , "fort.63"))
    np.testing.assert_allclose(d2, d1, rtol=rtol, atol=atol)
