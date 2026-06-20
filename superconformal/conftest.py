"""Put this folder on sys.path so `import sqcdkit` works when running
`pytest superconformal/tests/` from the repo root."""

import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
