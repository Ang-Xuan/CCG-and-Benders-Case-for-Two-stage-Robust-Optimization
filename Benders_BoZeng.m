clear all
%% Parameter
G=[-1,-1,-1,0,0,0,0,0,0;
    0,0,0,-1,-1,-1,0,0,0;
    0,0,0,0,0,0,-1,-1,-1;
    1,0,0,1,0,0,1,0,0;
    0,1,0,0,1,0,0,1,0;
    0,0,1,0,0,1,0,0,1];
h=[0;0;0;206;274;220];
E=[0,0,0,1,0,0;
    0,0,0,0,1,0;
    0,0,0,0,0,1;
    0,0,0,0,0,0;
    0,0,0,0,0,0;
    0,0,0,0,0,0];
M=[0,0,0;
    0,0,0;
    0,0,0;
    -40,0,0;
    0,-40,0;
    0,0,-40];
coe1 = [400,414,326];
coe2 = [18,25,20];
b = [22,33,24,33,23,30,20,25,27]' ;
%% Variable
y=binvar(3,1);
z=sdpvar(3,1);
x=sdpvar(9,1);
eta=sdpvar(1);
g=sdpvar(3,1);
pi = sdpvar(size(G,1),1);
v=binvar(size(G,1),1);
w=binvar(size(G,2),1);
%% Benders
LB=-inf;   UB=inf;   iter=1;   BigM=1e5;

MP_Obj = coe1*y +coe2*z+eta ;
MP_Cons=[ z<=800*y, z>=0, x>=0, eta>=b'*x, sum(z)>=772 ];
ops=sdpsettings('solver','cplex','verbose',0);

Uncertain_Cons=[ 0<=g<=1, sum(g)<=1.8, g(1)+g(2)<=1.2 ];

while UB-LB>=1e-5
    disp(['迭代第',num2str(iter),'次'])    
    optimize(MP_Cons,MP_Obj,ops);
    LB = max(LB, value(MP_Obj));                % LB

    SP_Obj = b'*x;
    SP_Cons = [Uncertain_Cons, 0<=x, G*x >= h-E*[value(y); value(z)]-M*g];
    SP_Cons = [SP_Cons, 0<=pi,  G'*pi<=b ];
    SP_Cons = [SP_Cons, G*x-h+E*[value(y); value(z)]+M*g <= BigM*(1-v), pi<=BigM*v];
    SP_Cons = [SP_Cons, b-G'*pi <= BigM*(1-w), x<=BigM*w];
    sol_SP = optimize(SP_Cons,-SP_Obj,ops);
    
    if sol_SP.problem==0                 % SP is solved
        UB=min(UB, coe1*value(y)+coe2*value(z)+value(SP_Obj));           % UB
        disp([' g =   ',num2str(value(g)')]);
        MP_Cons = [MP_Cons, eta>=( h-E*[y;z]-M*value(g) )'*value(pi)];
    else                                            % SP is unbounded
        MP_Cons = [MP_Cons];
    end
    
    iter=iter+1;
    display([' LB: ',num2str(LB), '    UB: ',num2str(UB),]);
end