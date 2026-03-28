function [c_TIT] = tit(state, eng)
    TIT     = state.TIT;
    TIT_max = eng.TIT_max;

    c_TIT   = (TIT - TIT_max) / TIT_max;
end