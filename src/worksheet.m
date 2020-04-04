%%
clear
clc
filename = "mesh/Case_A_param.input";

[M,P] = importInput(filename);

SP = feSolveParameters();
sol = feSolve(M,P,SP).solve()

%%
clear
clc
filename = "mesh/Case_B_param.input";

[M,P] = importInput(filename);

SP = feSolveParameters();
sol = feSolve(M,P,SP).solve()

%%
clear
clc
filename = "mesh/Case_C_param.input";

[M,P] = importInput(filename);

SP = feSolveParameters();
sol = feSolve(M,P,SP).solve()