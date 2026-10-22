use modena_core::{InputConfig,Simulation};
use std::path::Path;
fn main(){
 let input=include_str!("alloy.toml");
 let cfg=InputConfig::from_str(input).expect("config should pass strict validation");
 println!("InputConfig::from_str accepted valid positive concentrations totaling1");
 match std::panic::catch_unwind(||Simulation::from_config_without_outputs(&cfg,Path::new("."))){
  Ok(Ok(s))=>println!("constructed {} atoms",s.system().n()),
  Ok(Err(e))=>println!("returned Result error: {e}"),
  Err(e)=>{
   let msg=e.downcast_ref::<&str>().copied().or_else(||e.downcast_ref::<String>().map(String::as_str)).unwrap_or("unknown panic");
   println!("PANIC instead of Result: {msg}");
  }
 }
}
