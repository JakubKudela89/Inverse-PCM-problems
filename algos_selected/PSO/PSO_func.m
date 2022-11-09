function [res]= PSO_func(Dimension,Particle_Number,fhd,Max_Gen,VRmin,VRmax)
%[gbest,gbestval,fitcount]= PSO_func('f8',3500,200000,30,30,-5.12,5.12)
rand('state',sum(100*clock));
me=Max_Gen;
ps=Particle_Number;
D=Dimension;
cc=[2 2];   %acceleration constants
iwt=0.9-(1:me).*(0.5./me);
% iwt=0.5.*ones(1,me);
if length(VRmin)==1
    VRmin=repmat(VRmin,1,D);
    VRmax=repmat(VRmax,1,D);
end
mv=0.5*(VRmax-VRmin);
VRmin=repmat(VRmin,ps,1);
VRmax=repmat(VRmax,ps,1);
Vmin=repmat(-mv,ps,1);
Vmax=-Vmin;
pos=VRmin+(VRmax-VRmin).*rand(ps,D);

% for n=0:9
%     RecordFEsFactor(n+1) = round(D^((n/10)-0.9)*Max_Gen*ps);
% end
nr_points_to_save = 10;
RecordFEsFactor = round(logspace(log(2*ps)/log(10),log(me*ps)/log(10),nr_points_to_save));
Run_RecordFEsFactor=RecordFEsFactor;
current_eval = 0; 

e = zeros(1,ps);
e = [fhd(pos(1:ps,:))]';

run_funcvals = [];
run_funccalls = [];

%e=feval(fhd,pos',varargin{:});

fitcount=ps;
vel=Vmin+2.*Vmax.*rand(ps,D);%initialize the velocity of the particles
pbest=pos;
pbestval=e; %initialize the pbest and the pbest's fitness value
[gbestval,gbestid]=min(pbestval);
gbest=pbest(gbestid,:);%initialize the gbest and the gbest's fitness value
gbestrep=repmat(gbest,ps,1);
current_eval = current_eval + Particle_Number;

run_funcvals = [run_funcvals; gbestval];
run_funccalls = [run_funccalls; current_eval];
fprintf('fes: %u, val: %5.3e\n', current_eval, gbestval)
            
for i=2:me
        aa=cc(1).*rand(ps,D).*(pbest-pos)+cc(2).*rand(ps,D).*(gbestrep-pos);
        vel=iwt(i).*vel+aa;
        vel=(vel>Vmax).*Vmax+(vel<=Vmax).*vel;
        vel=(vel<Vmin).*Vmin+(vel>=Vmin).*vel;
        pos=pos+vel;
        pos=((pos>=VRmin)&(pos<=VRmax)).*pos...
            +(pos<VRmin).*(VRmin+0.25.*(VRmax-VRmin).*rand(ps,D))+(pos>VRmax).*(VRmax-0.25.*(VRmax-VRmin).*rand(ps,D));
        %e=feval(fhd,pos',varargin{:});
        e = zeros(1,ps);
        e = [fhd(pos(1:ps,:))]';

        fitcount=fitcount+ps;
        tmp=(pbestval<e);
        temp=repmat(tmp',1,D);
        pbest=temp.*pbest+(1-temp).*pos;
        pbestval=tmp.*pbestval+(1-tmp).*e;%update the pbest
        [gbestval,tmp]=min(pbestval);
        gbest=pbest(tmp,:);
        gbestrep=repmat(gbest,ps,1);%update the gbest
        current_eval = current_eval + Particle_Number;
        run_funcvals = [run_funcvals; gbestval];
        run_funccalls = [run_funccalls; current_eval];
        fprintf('fes: %u, val: %5.3e\n', current_eval, gbestval)
        
        if gbestval<1e-8
            run_funcvals(end+1:36) = 0;
           break 
        end
end
    
if gbestval < 1e-8
   gbestval = 0; 
end

RecordFEsFactor(RecordFEsFactor<1e-8) = 0;

run_funcvals(end) = gbestval;
res.solution = gbest;
res.bestval = gbestval;
res.progress_iters = run_funccalls;
res.progress_values = run_funcvals;

end


