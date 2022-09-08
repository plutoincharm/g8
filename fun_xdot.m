function xdot = fun_xdot(x,u,dt)

global A B Nx Nu pert MI L m  nx ny tx ty g ra lam vars alp alpval indic kc lamall xdata lamx lamy

xdot = [x(10:18); fun_qddotsspsinglecontact(x,u,dt)];
%xdot = fun_qddotssp11(x,u,dt);
