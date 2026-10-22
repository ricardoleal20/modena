use modena_core::{System,SimBox,Vec3,ForceModel,Virial,NeighborList,SpmeForceModel,EwaldForceModel,InputConfig,Simulation};
use std::path::Path;
fn eval(m:&dyn ForceModel,s:&System)->(f64,Vec3){
 let mut x=s.clone(); x.forces.fill(0.0); let mut v=Virial::default();let mut n=NeighborList::with_skin(4.2,0.3);
 let u=m.compute(&mut x,&mut v,&mut n).unwrap();
 let mut f=Vec3::zero();for i in 0..x.n(){f+=x.force(i);} (u,f)
}
fn momentum(s:&System)->Vec3{let mut p=Vec3::zero();for i in 0..s.n(){p+=s.vel(i)*s.masses[i];}p}
fn main(){
 let mut s=System::new(SimBox::cubic(9.0));s.resize(3);s.charges=vec![1.0,-0.7,-0.3];
 for(i,p)in [Vec3::new(1.1,2.3,3.4),Vec3::new(2.7,3.2,3.8),Vec3::new(6.2,5.8,7.0)].into_iter().enumerate(){s.set_pos(i,p);}
 for k in [3.2,4.8,6.4]{
 let m=SpmeForceModel::new(0.8,4.2,k).unwrap();let(u,f)=eval(&m,&s);let mut hi=s.clone();let mut lo=s.clone();let h=1e-5;
 for i in 0..s.n(){hi.set_pos(i,s.pos(i)+Vec3::new(h,0.0,0.0));lo.set_pos(i,s.pos(i)-Vec3::new(h,0.0,0.0));}
 println!("spme kcut={k} U={u:.16e} netF={:?} minus_dU_dx={:.16e}",f.0,-(eval(&m,&hi).0-eval(&m,&lo).0)/(2.0*h));
 }
 let direct=EwaldForceModel::new(0.8,4.2,6.4).unwrap();println!("ewald netF={:?}",eval(&direct,&s).1.0);
 let base=Path::new(".superpowers/sdd/2026-10-03-scientific-hardening/final-probe");
 for method in ["spme","ewald"]{
 let raw=format!(r#"units="reduced"
[structure]
type="xyz"
file="charges.xyz"
box_size=[9.0,9.0,9.0]
[model]
type="lennard_jones"
sigma=1.0
epsilon=0.0
cutoff=4.2
[integrator]
type="velocity_verlet"
dt=0.001
[electrostatics]
method="{method}"
kappa=0.8
real_cutoff=4.2
kspace_cutoff=4.8
[run]
steps=100
initial_temperature=1.0
seed=51
log_every=100
"#);
 let cfg=InputConfig::from_str(&raw).unwrap();let mut sim=Simulation::from_config_without_outputs(&cfg,base).unwrap();
 println!("runner {method} initial metadata={:?} P={:?}",sim.metadata(),momentum(sim.system()).0);
 sim.run().unwrap();println!("runner {method} final metadata={:?} P={:?}",sim.metadata(),momentum(sim.system()).0);
 }
}
