import os
from pathlib import Path

__all__ = ['DRIVE', 'DROPBOX', 'drive_path']

def drive_path():
	# Root holding research_files/ and documents/: the Google Drive "My Drive" folder.
	# Override with the DRIVE_ROOT env var (e.g. on the cluster, where Google Drive is not mounted).
	if 'DRIVE_ROOT' in os.environ:
		return Path(os.environ['DRIVE_ROOT'])
	# Google Drive for desktop mounts each account at ~/Library/CloudStorage/GoogleDrive-<email>/My Drive;
	# pick the account with GDRIVE_ACCOUNT (its email) if there is more than one.
	account = os.environ.get('GDRIVE_ACCOUNT', '*')
	drives = sorted((Path.home() / 'Library' / 'CloudStorage').glob(f'GoogleDrive-{account}/My Drive'))
	return drives[0] if drives else Path.home() / 'Google Drive' / 'My Drive'

DRIVE = drive_path()
DROPBOX = DRIVE   # deprecated alias, kept for scripts outside this repository
