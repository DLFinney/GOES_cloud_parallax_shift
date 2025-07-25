#!/usr/bin/env python3
"""
Example usage of GOES cloud parallax shift tools

This script demonstrates the basic workflow for downloading GOES data
and applying parallax correction.
"""

import os
import sys
import subprocess

def run_command(cmd, description):
    """Run a command and handle errors"""
    print(f"\n{description}")
    print(f"Running: {cmd}")
    
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print("SUCCESS")
        if result.stdout:
            print(f"Output: {result.stdout}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: {e}")
        if e.stderr:
            print(f"Error details: {e.stderr}")
        return False

def main():
    """Main example workflow"""
    print("GOES Cloud Parallax Shift - Example Usage")
    print("=" * 50)
    
    # Check if we're in the right directory
    if not os.path.exists("download_goes_subregion_regrid.py"):
        print("Error: Please run this script from the repository root directory")
        sys.exit(1)
    
    print("\nThis example demonstrates the basic workflow:")
    print("1. Download GOES satellite data")
    print("2. Apply parallax correction")
    print("\nNote: This is a demonstration only. Actual data download requires")
    print("internet connection and may take significant time.")
    
    # Example 1: Download a single product for one month
    cmd1 = "python download_goes_subregion_regrid.py ABI-L2-ACHAC 7"
    run_command(cmd1, "Step 1: Download cloud top height data for July")
    
    # Example 2: Apply parallax correction
    cmd2 = "python parallax_latlons_goes_cloud.py"
    run_command(cmd2, "Step 2: Apply parallax correction")
    
    print("\n" + "=" * 50)
    print("Example complete!")
    print("\nFor actual usage:")
    print("1. Review and modify the configuration sections in each script")
    print("2. Set appropriate date ranges, regions, and output paths")
    print("3. Ensure you have sufficient disk space for satellite data")
    print("4. Check that all dependencies are properly installed")

if __name__ == "__main__":
    main()