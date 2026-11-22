x1 = 0.0; x2 = 1.0;
y1 = 0.0; y2 = 12.0;

x1_px = getResult("X", 0); x2_px = getResult("X", 1);
y1_px = getResult("Y", 2); y2_px = getResult("Y", 3);

n = nResults;
for (i = 4; i < n; i++) {
    t_us = x1 + (getResult("X", i) - x1_px) / (x2_px - x1_px) * (x2 - x1);
    R_um = y1 + (getResult("Y", i) - y1_px) / (y2_px - y1_px) * (y2 - y1);
    setResult("dt_us", i, t_us);
    setResult("R_um", i, R_um);
}
updateResults();