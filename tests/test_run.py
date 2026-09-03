#!/usr/bin/env python3
import pytest
from util import run_serial, run_parallel, run_dagswem
import os


def test_quarter_annular(binpath, test_dir):
    run_serial(binpath, str(test_dir("quarter_annular")), 0.05, 1e-7)

def test_quarter_annular_parallel(binpath, test_dir):
    run_parallel(binpath, str(test_dir("quarter_annular")), 0.05, 1e-7)

def test_quarter_annular_dagswem(binpath, test_dir):
    run_dagswem(binpath, str(test_dir("quarter_annular")), 0.05, 1e-7, num_ranks=2)

def test_performance_quarter_annular(binpath, test_dir, mpi_aps):
    run_parallel(binpath, str(test_dir("quarter_annular")), 0.05, 1e-7, num_ranks=4)

def test_wetdry(binpath, test_dir):
    run_serial(binpath, str(test_dir("wetdry")), 0.01, 1e-7)

def test_wetdry_parallel(binpath, test_dir):
    run_parallel(binpath, str(test_dir("wetdry")), 0.01, 1e-7)

def test_mass_conservation(binpath, test_dir):
    run_serial(binpath, str(test_dir("mass_conservation")), 1e-3, 1e-7)

