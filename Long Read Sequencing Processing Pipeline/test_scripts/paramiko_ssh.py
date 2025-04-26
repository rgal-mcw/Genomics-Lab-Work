import paramiko

def create_ssh_client(server, port, user, key_path):
    """Create and return an SSH client, connected to the server."""
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    key = paramiko.RSAKey.from_private_key_file(key_path)
    client.connect(server, port, user, pkey=key)
    return client

def run_command_via_ssh(client, command):
    """Run a command on the server, using the SSH client."""
    stdin, stdout, stderr = client.exec_command(command)
    output = stdout.read()
    error = stderr.read()
    
    stdin.close()
    stdout.close()
    stderr.close()

    if error:
        return "Error: " + error.decode()
    return output.decode()


def main():
    server = 'maple.phys.mcw.edu'
    port = 22
    username = 'rgallagher'
    key = '/Users/ry28926/.ssh/id_rsa'

    ssh_client = create_ssh_client(server, port, username, key)

    print("Running 'ls -l' on the server:")
    output = run_command_via_ssh(ssh_client, 'ls -l')
    print("Output:")
    print(output)

    ssh_client.close()

if __name__ == '__main__':
    main()

