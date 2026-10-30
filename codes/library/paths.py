import json
import os
from pathlib import Path

__all__ = ['DROPBOX', 'dropbox_path']

def dropbox_path(account=None):
	# Pick the Dropbox account with DROPBOX_ACCOUNT ('personal' or 'business');
	# defaults to 'personal'.
	account = account or os.environ.get('DROPBOX_ACCOUNT', 'business')
	# Dropbox records its folder location in ~/.dropbox/info.json on every install;
	# fall back to ~/Dropbox (folder or symlink) if that file is missing.
	try:
		with open(Path.home() / '.dropbox' / 'info.json') as f:
			return Path(json.load(f)[account]['path'])
	except (OSError, KeyError, ValueError):
		return Path.home() / 'Dropbox'

DROPBOX = dropbox_path()
