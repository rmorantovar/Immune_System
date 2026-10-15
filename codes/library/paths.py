import json
from pathlib import Path

__all__ = ['DROPBOX', 'dropbox_path']

def dropbox_path():
	# Dropbox records its folder location in ~/.dropbox/info.json on every install;
	# fall back to ~/Dropbox (folder or symlink) if that file is missing.
	try:
		with open(Path.home() / '.dropbox' / 'info.json') as f:
			return Path(json.load(f)['personal']['path'])
	except (OSError, KeyError, ValueError):
		return Path.home() / 'Dropbox'

DROPBOX = dropbox_path()
