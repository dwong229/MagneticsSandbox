function statedot = compute_magnet_deriv_4coil(t,state,magnet,coil,params)

x = state(1:2);
v = state(3:4);

% Compute F :
a = computeF_4coil(coil,magnet,params,x,v)/magnet.mass;

% pack output
statedot = [v;a];







