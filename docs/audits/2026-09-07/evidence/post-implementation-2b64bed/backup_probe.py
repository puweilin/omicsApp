"""Run from the reviewed repository root; uses only a temporary docker stub."""
from pathlib import Path
import datetime
import os
import subprocess
import tempfile

with tempfile.TemporaryDirectory(prefix="omicsapp-backup-review-") as temp:
    root = Path(temp)
    (root / "bin").mkdir()
    (root / "store").mkdir()
    docker = root / "bin" / "docker"
    docker.write_text(
        '#!/bin/bash\n'
        'printf "PostgreSQL database dump\\n"\n'
        '/usr/bin/awk \'BEGIN {for(i=0;i<200000;i++) '
        'print "INSERT INTO users VALUES (123, 456, 789, 123456789);"}\'\n'
    )
    docker.chmod(0o755)
    day = datetime.date.today().isoweekday()
    destination = root / "store" / f"keycloak-{day}.sql.gz"
    previous = b"previous-backup-preserve-me"
    destination.write_bytes(previous)
    environment = os.environ.copy()
    environment["PATH"] = str(root / "bin") + ":" + environment["PATH"]
    environment["BACKUP_DIR"] = str(root / "store")
    result = subprocess.run(
        ["bash", "deploy/scripts/backup-keycloak.sh"],
        env=environment, capture_output=True, text=True,
    )
    print(f"LARGE_VALID_DUMP_EXIT={result.returncode}")
    print(f"PREVIOUS_BACKUP_PRESERVED={destination.read_bytes() == previous}")
    print(f"STDERR={result.stderr}")
