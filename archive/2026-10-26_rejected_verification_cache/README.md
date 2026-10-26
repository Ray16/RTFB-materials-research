# Rejected CPU verification-cache experiment

This is preserved failed work, excluded from active workflow code and test discovery.
The cache failed two corruption tests: fast same-size writes could leave timestamp
metadata unchanged, including a write that restored mtime. It was briefly enabled
before those test results were checked. The CPU controller was stopped and the
unmodified strict coordinator restored. Scientific worker source checks and
numerical code were never changed. No source data were altered by these tests.
The failure motivated retaining complete checks; this archive makes no performance
or integrity qualification claim. Original paths and hashes are in manifest.json.
