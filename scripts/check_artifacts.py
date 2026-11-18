"""Validate the exact certificate tape and the reviewed manuscript artifact."""

import base64
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_hashes(records):
    for name, expected in records.items():
        if sha256(ROOT / name) != expected:
            raise ValueError(f"SHA-256 mismatch: {name}")


def main():
    manifest = json.loads((ROOT / "certificate/manifest.json").read_text(encoding="utf-8"))
    verify_hashes(manifest["verification_payloads"])
    tape = (ROOT / "certificate/three-star-tape.txt").read_text(encoding="ascii")
    lines = tape.splitlines()
    if lines[0] != "CERTIFICATE_BEGIN" or lines[-1] != "CERTIFICATE_END":
        raise ValueError("Certificate markers are missing.")
    raw = base64.b64decode("".join(lines[1:-1]), validate=True)
    if len(raw) != manifest["decoded_bytes"] or hashlib.sha256(raw).hexdigest() != manifest["decoded_sha256"]:
        raise ValueError("Decoded certificate does not match the manifest.")
    paper = json.loads((ROOT / "manuscript/manifest.json").read_text(encoding="utf-8"))
    verify_hashes(paper["sha256"])
    print(f"Verified {len(raw):,} certificate bytes, both checker hashes, and the reviewed manuscript files.")


if __name__ == "__main__":
    main()
