% File: hammerDeltaL.m
function dL = hammerDeltaL(eeff, h, W)
dL = 0.412*h .* ((eeff+0.3).*(W/h+0.264)) ./ ((eeff-0.258).*(W/h+0.8));
end