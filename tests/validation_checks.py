# validation_checks.py
import sys
from pathlib import Path

def validate_paths(database_path, file_path):
    """
    The validate function where it checks for valid paths.
    """
   
    database_path = Path(database_path)
    file_path = Path(file_path)
    
    paths_to_check = {
        'Database Path': database_path,
        'File Path': file_path
    }

    for name, path in paths_to_check.items():
        if not path.exists():
            print(f"Error: {name} not found.", file=sys.stderr)
            print(f"       Location checked: {path}", file=sys.stderr)
            sys.exit(1)

    # print("Success! Both database path and file path exist.")
    # print(f"Database: {database_path}")
    # print(f"File: {file_path}")
    return