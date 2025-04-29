# Title: ONT -> Maple Data Transfer
# Author: Ryan Gallagher (Broeckel Lab, 2024)
# 
# Descriptions:
#   This script creates an ssh connection from Maple to the PromethION, then
#   runs various file management commands before transfering the data to Maple.
#
#   The script is designed to be run from the Maple server, and will prompt then
#   user for the necessary information about file management within Maple.
#
#   This script:
#       1. Asks for the directory on the prom where the file is stores
#       2. Creates an ssh connection to the prompt
#       3. Runs various file management commands
#       4. Closes the ssh connection
#       5. Transfers the data to Maple (from the maple session)
#
#
#   IMPORTANT: THIS SCRIPT REQUIRES AN SSH KEY TO BE INSTALLED ON THE PROMETHION
#
#   CW Maple Transfer Revisit (04.2025):
#       * PROM now has a new IP address
#           -> Changed from 141.106.247.152 to 10.41.1.34
#           -> The port should not have changed. I think this is a prom thing.
#       * Updated create_ssh_client() function to look for rsa key instead of inputting password
#       * Included sys arguments for better error handling

import os
import paramiko
import subprocess
import argparse
import sys
import logging

def setup_logging(log_file_path, log_level=logging.INFO): # <-- Default log level
    """Configures logging for the script."""
    # Use basicConfig to configure the root logger
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file_path, mode='a'), # Append mode
            logging.StreamHandler(sys.stdout) # Log to console
        ]
    )
    # Return a logger specific to this module
    return logging.getLogger(__name__)


def create_ssh_client(server, port, user, logger):
    """
    Create and return an SSH client, connected to the server using an RSA key.

    Assumes the private key is located at ~/.ssh/id_rsa.
    """
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    # Specify the path to the RSA private key file
    key_path = os.path.expanduser('~/.ssh/id_rsa') # Using standard RSA key path
    logger.info(f"Attempting SSH connection to {user}@{server}:{port} using key: {key_path}")
    try:
        # Use the key_filename parameter for key-based authentication
        # No password parameter needed when using key authentication
        client.connect(server, port, user, key_filename=key_path)
        logger.info("SSH connection established successfully.")
        return client
    except paramiko.AuthenticationException:
        logger.error("Authentication failed duing SSH connection.")
        print("\nAuthentication failed.")
        print("Please check:")
        print(f"  - If the public key corresponding to:\n    {key_path}")
        print("    is in the authorized_keys file on the remote server for user '{user}'.")
        print("  - If your local private key has the correct permissions (usually 600).")
        # Note: If you are using a passphrase, paramiko might prompt for it interactively
        return None
    except paramiko.SSHException as e:
        logger.error(f"SSH connection failed: {e}")
        print(f"\nSSH connection failed: {e}")
        print("Please check:")
        print(f"  - If the server '{server}' is reachable.")
        print(f"  - If port '{port}' is open on the server.")
        return None
    except FileNotFoundError:
         logger.error(f"SSH key not found at {key_path}.")
         print(f"\nError: SSH key not found at {key_path}.")
         print("Please ensure the private key file exists at this path.")
         return None
    except Exception as e:
        logger.exception("An unexpected error occured during SSH connection.")
        print(f"\nAn unexpected error occurred during SSH connection: {e}")
        return None


def check_directory_exists(ssh_client, directory, logger):
    """Check if a directory exists on the remote server via SSH."""
    # Use a command that returns 0 for success, non-zero for failure
    command = f'test -d "{directory}"'
    logger.info(f"Checking remote directory: {command}")
    stdin, stdout, stderr = ssh_client.exec_command(command)

    # Read output and errors to prevent hanging
    stdout.read()
    stderr_output = stderr.read().decode().strip()


    exit_status = stdout.channel.recv_exit_status() # Get the exit status of the remote command

    if exit_status != 0:
        # Directory does not exist or other error occurred
        if stderr_output:
             logger.warning(f"Remote test command stderr: {stderr_output}")
        return False
    return True # Directory exists (exit status was 0)

def run_command_via_ssh(client, command, logger):
    """Run a command on the server, using the SSH client, and return output/errors."""
    logger.info(f"Executing remote command: {command}")
    try:
        stdin, stdout, stderr = client.exec_command(command)

        # Read output and error streams completely
        output = stdout.read().decode('utf-8', errors='ignore').strip()
        error = stderr.read().decode('utf-8', errors='ignore').strip()

        # Wait for the command to finish and get exit status
        exit_status = stdout.channel.recv_exit_status()

        # Check if the remote command failed (non-zero exit status)
        if exit_status != 0:
            logger.warning(f"Remote command failed with exit status {exit_status}")
            if error:
                logger.error("Remote command Error:\n" + error)
            # Indicate failure by returning an error string or raising exception
            return f"ERROR: Remote command failed (exit status {exit_status}). Stderr:\n{error}\nStdout:\n{output}"

        # Print stderr even if command succeeded (might contain warnings)
        if error:
             logger.error("Remote command Stderr (warnings/non-fatal issues):\n" + error)

        return output

    except Exception as e:
        logger.exception(f"An error occurred while executing remote command via SSH: {e}")
        # Indicate failure
        return f"ERROR: Error executing command via SSH: {e}"



def run_rsync(source, destination, logger, options=['-a', '-e', 'ssh -i ~/.ssh/id_rsa']):
    """
    Run an rsync command with specified source, destination, and options.

    Default options include archive mode (-a) and using the RSA SSH key (-e 'ssh -i ~/.ssh/id_rsa').
    """
    rsync_command = ["rsync"] + options + [source, destination]
    # Log the command being run (NEW)
    logger.info("Running rsync command: " + " ".join(rsync_command))

    try: # NEW: Added try block for error handling
        # Use check=True to raise a CalledProcessError on non-zero exit codes (NEW)
        result = subprocess.run(rsync_command, capture_output=True, text=True, check=True)
        logger.info("Rsync completed successfully.")
        # Keep original behavior, only print stderr on error (Output printing removed)
        # if result.stdout:
        #     print("Rsync Stdout:\n", result.stdout)

    except subprocess.CalledProcessError as e: # NEW: Specific handling for command failure
        logger.error("\nError during rsync.")
        logger.error("Rsync command failed with exit status:", e.returncode)
        if e.stderr:
            logger.error("Rsync Stderr:\n", e.stderr) # Print stderr from the exception
        # if e.stdout: # Don't print stdout on error unless specifically needed
        #      print("Rsync Stdout:\n", e.stdout)
        # Exit the script with a non-zero status on rsync failure (NEW)
        sys.exit(f"Rsync failed with exit code {e.returncode}")
    except FileNotFoundError: # NEW: Specific handling if rsync command is not found
         logger.error("Error: The rsync command was not found.")
         print(f"\nError: The rsync command was not found.")
         print("Please ensure 'rsync' is installed and in your system's PATH.")
         sys.exit("Rsync command not found") # Exit the script
    except Exception as e: # NEW: Catch any other unexpected errors
         logger.exception("An unexpected error occured during rsync:")
         print(f"\nAn unexpected error occurred during rsync: {e}")
         sys.exit(f"An unexpected error occurred: {e}") # Exit the script


def main():
    # Define the ssh target
    server = '10.41.1.34'
    port = 22
    username = 'prom'

    # Take inputs from the main pipeline.py file
    parser = argparse.ArgumentParser(description='Process inputs')
    parser.add_argument('sample_name', type=str, help='Enter the sample name')
    parser.add_argument('data_dir', type=str, help='Enter the directory where the data is stored on the PromethION')
    parser.add_argument('to_dir', type=str, help='Enter the MAPLE directory where this should be stored')
    args = parser.parse_args()

    if not args.data_dir.endswith('/'):
        args.data_dir += '/'

    if not args.to_dir.endswith('/'):
        args.to_dir +='/'

    prom_name = f'{username}@{server}'

    log_dir = os.path.join(args.to_dir, 'logs')
    os.makedirs(log_dir, exist_ok=True)
    log_file = os.path.join(log_dir, 'file_transfer.log')
    # Assign the return value of setup_logging to the 'logger' variable
    logger = setup_logging(log_file, log_level=logging.INFO)

    # Write to log
    logger.info("The sample name is: " + args.sample_name)
    logger.info("The PROM data directory is: " + args.data_dir)
    logger.info("The MAPLE data directory is: " + args.to_dir)

    # Start ssh using credentials
    # Called without password since we're using id_rsa
    ssh_client = create_ssh_client(server, port, username, logger)
    
    if ssh_client is None:
        logger.error("Failed to establish SSH connection. Exiting.")
        sys.exit() # Using original exit()

    if check_directory_exists(ssh_client, args.data_dir, logger):
        logger.info("The PROM directory exists.")
    else:
        logger.error("The PROM directory does not exist. Exiting script.")
        ssh_client.close()
        sys.exit(1)
    # File management via ssh
    command = (

        f"cd {args.data_dir}",
        '[ -n "$(ls *.txt 2>/dev/null)" ] && mv *.txt other_reports || echo "No .txt files to move"',
        '[ -n "$(ls *.csv 2>/dev/null)" ] && mv *.csv other_reports || echo "No .csv files to move"',
        '[ -n "$(ls *.html 2>/dev/null)" ] && mv *.html other_reports || echo "No .html files to move"',
        '[ -n "$(ls *.json 2>/dev/null)" ] && mv *.json other_reports || echo "No .json files to move"',
        '[ -n "$(ls *.md 2>/dev/null)" ] && mv *.md other_reports || echo "No .md files to move"',
        '[ -n "$(ls *.tsv 2>/dev/null)" ] && mv *.tsv other_reports || echo "No .tsv files to move"',
        f"cat ./fastq_pass/*.fastq.gz > {args.sample_name}.fastq.gz",
        
        # Zipping pod_5 will not reduce the size of the directory. (Already in binary)
        #'tar -zcvf pod5_pass.tar.gz pod5_pass'
    )

    command = " && ".join(command)

    # Log and run command
    logger.info(f"Running command: {command}")
    output = run_command_via_ssh(ssh_client, command, logger)
    logger.info("Output:")
    logger.info(f"{output}")
    
    # Close ssh
    logger.info("Closing SSH connection.")
    ssh_client.close()

    # Begin rsync from prom to maple
    logger.info("Starting rsync data transfer to Maple...")
    destination1 = f"{prom_name}:{args.data_dir}{args.sample_name}.fastq.gz"
    to_dir1 = f"{args.to_dir}raw_data/"
    if not os.path.exists(to_dir1):
        os.makedirs(to_dir1)

    run_rsync(destination1, to_dir1, logger)
    logger.info(f"Transferred {args.sample_name}.fastq.gz to {to_dir1}")

    destination2 = f"{prom_name}:{args.data_dir}other_reports"
    
    run_rsync(destination2, to_dir1, logger)
    logger.info(f"Transferred other_reports to {to_dir1}")

    destination3 = f"{prom_name}:{args.data_dir}pod5"

    logger.info("Transferring pod5 directory")
    run_rsync(destination3, to_dir1, logger)
    logger.info(f"Transferred pod5/ to {to_dir1}")
    
    #print("pod5_pass.tar.gz NOT COPIED OR GZIP'D for storage reason. Edit the script to copy if needed.")


if __name__ == '__main__':
    main()

