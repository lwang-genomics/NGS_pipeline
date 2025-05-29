import os
import sys
import time
import subprocess
import re

def get_sample_name(filename):
    base = os.path.basename(filename)
    match = re.match(r"(.+?)\.R[12]\.f(?:ast)?q(?:\.gz)?$", base)
    if not match:
        raise ValueError(f"Invalid filename format: {base}. Expected format: sample.R1.fq or sample.R1.fastq(.gz)")
    return match.group(1)

def start_logging(sample_name):
    log_filename = f"{sample_name}.log"
    log_file = open(log_filename, "w")

    start_time = time.strftime("%Y-%m-%d %H:%M:%S")
    border = "=" * 100

    log_file.write(f"\n{border}\n")
    log_file.write(f"{'[PIPELINE STARTED]'.center(100)}\n")
    log_file.write(f"{border}\n")
    log_file.write(f"{'Start Time:'.ljust(20)} {start_time}\n")
    log_file.write(f"{'Executed Command:'.ljust(20)} {' '.join(sys.argv)}\n")
    log_file.write(f"{border}\n\n")
    log_file.flush()

    return log_file

def log_stage(stage_name, command, log_file):
    border = "-" * 100
    log_file.write(f"\n{border}\n")
    log_file.write(f"{f'STAGE: {stage_name}':^100}\n")
    log_file.write(f"{border}\n")
    log_file.write(f"{'[COMMAND]'.ljust(15)} {command}\n")
    log_file.flush()

    try:
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        if result.stdout:
            log_file.write(f"\n{'[STDOUT]'.ljust(15)}\n{result.stdout}\n")
        if result.stderr:
            log_file.write(f"\n{'[STDERR]'.ljust(15)}\n{result.stderr}\n")
        if result.returncode != 0:
            log_file.write(f"\n{'[ERROR]'.ljust(15)} Stage '{stage_name}' failed with return code {result.returncode}\n")
            log_file.write(f"{border}\n")
            log_file.flush()
            log_file.close()
            sys.exit(result.returncode)
    except Exception as e:
        log_file.write(f"\n{'[EXCEPTION]'.ljust(15)} {str(e)}\n")
        log_file.write(f"{border}\n")
        log_file.flush()
        log_file.close()
        sys.exit(1)

    log_file.write(f"\n{'[STATUS]'.ljust(15)} Stage '{stage_name}' completed successfully\n")
    log_file.write(f"{border}\n")
    log_file.flush()

def end_logging(log_file):
    end_time = time.strftime("%Y-%m-%d %H:%M:%S")
    border = "=" * 100
    log_file.write(f"\n{border}\n")
    log_file.write(f"{'[PIPELINE COMPLETED]'.center(100)}\n")
    log_file.write(f"{border}\n")
    log_file.write(f"{'End Time:'.ljust(20)} {end_time}\n")
    log_file.write(f"{border}\n")
    log_file.close()