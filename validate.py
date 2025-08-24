#!/usr/bin/env python3
"""
Validation script for GOES cloud parallax shift tools

This script performs basic validation to ensure the tools are working correctly.
"""

import os
import sys
import importlib
import subprocess

def test_import(module_name, friendly_name=None):
    """Test if a module can be imported"""
    if friendly_name is None:
        friendly_name = module_name
    
    try:
        importlib.import_module(module_name)
        print(f"✓ {friendly_name} - OK")
        return True
    except ImportError as e:
        print(f"✗ {friendly_name} - MISSING ({e})")
        return False

def test_script_syntax(script_path):
    """Test if a Python script has valid syntax"""
    try:
        result = subprocess.run([sys.executable, '-m', 'py_compile', script_path], 
                              capture_output=True, text=True)
        if result.returncode == 0:
            print(f"✓ {script_path} - Syntax OK")
            return True
        else:
            print(f"✗ {script_path} - Syntax Error: {result.stderr}")
            return False
    except Exception as e:
        print(f"✗ {script_path} - Error: {e}")
        return False

def main():
    """Run validation tests"""
    print("GOES Cloud Parallax Shift - Validation")
    print("=" * 40)
    
    # Check we're in the right directory
    if not os.path.exists("download_goes_subregion_regrid.py"):
        print("Error: Please run this script from the repository root directory")
        sys.exit(1)
    
    print("\n1. Testing core dependencies:")
    all_deps_ok = True
    
    # Test core Python modules
    deps = [
        ('numpy', 'NumPy'),
        ('xarray', 'xarray'),
        ('matplotlib', 'Matplotlib'),
        ('pandas', 'Pandas')
    ]
    
    for module, name in deps:
        if not test_import(module, name):
            all_deps_ok = False
    
    print("\n2. Testing optional dependencies:")
    # Test optional modules (these may fail and that's OK)
    optional_deps = [
        ('goes2go', 'goes2go (satellite data access)'),
        ('xesmf', 'xesmf (regridding)')
    ]
    
    for module, name in optional_deps:
        test_import(module, name)
    
    print("\n3. Testing script syntax:")
    scripts = [
        'download_goes_subregion_regrid.py',
        'parallax_latlons_goes_cloud.py',
        'example_usage.py',
        'validate.py'
    ]
    
    all_syntax_ok = True
    for script in scripts:
        if os.path.exists(script):
            if not test_script_syntax(script):
                all_syntax_ok = False
        else:
            print(f"✗ {script} - File not found")
            all_syntax_ok = False
    
    print("\n4. Testing file structure:")
    required_files = [
        'README.md',
        'requirements.txt',
        '.gitignore',
        'example_input_files/README',
        'example_output_files/README'
    ]
    
    structure_ok = True
    for file in required_files:
        if os.path.exists(file):
            print(f"✓ {file} - Present")
        else:
            print(f"✗ {file} - Missing")
            structure_ok = False
    
    print("\n" + "=" * 40)
    print("VALIDATION SUMMARY:")
    
    if all_deps_ok:
        print("✓ Core dependencies: OK")
    else:
        print("✗ Core dependencies: Some missing")
    
    if all_syntax_ok:
        print("✓ Script syntax: OK")
    else:
        print("✗ Script syntax: Errors found")
    
    if structure_ok:
        print("✓ File structure: OK")
    else:
        print("✗ File structure: Missing files")
    
    if all_deps_ok and all_syntax_ok and structure_ok:
        print("\n🎉 All validation tests passed!")
        return 0
    else:
        print("\n⚠️  Some validation tests failed. See details above.")
        return 1

if __name__ == "__main__":
    sys.exit(main())