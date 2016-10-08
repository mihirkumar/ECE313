function [out] = threshold_func(in, a, b)
    out = ((in < a) | (in > b));
