"""Run the full Step1 -> Step2 -> Step3 pipeline on a few examples."""

from conformalmanifold import run
from conformalmanifold.groups import cyclic

# Step 1 can take a library name ...
run("A4 = Delta(12)")
run("Delta(27)")
run("A5 = Sigma(60)")

# ... or a group you build yourself.
run(cyclic(10, (2, 3, 5)))

# Programmatic access to just the answer:
out = run("Delta(27)", verbose=False)
print("\nDelta(27): dim_C M_conf =", out["dim_conf"])
