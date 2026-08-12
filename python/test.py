from reproduce_figures import *
import numpy as np
lam = mg_to_lambda(1.6)
w, c, ns, mi, st, nso, Nm = run_dlcq(3, 1, 1, 21, lam)
x, qq, qa = extract_structure_function(0, c, nso, mi, st.mstate, Nm, 21, 3, 1)
print(f"Baryon sum rule: {np.sum(qq) * 2.0/21:.4f}  (should be 3.0)")
print(f"Peak: {np.max(qq):.2f} at x={x[np.argmax(qq)]:.3f}")
