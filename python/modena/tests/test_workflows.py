import math
import os
from pathlib import Path
import sys
import tempfile
import tomllib
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import modena.config as c
from modena import check,run
BINARY=os.environ.get("MODENA_BINARY",str(Path(__file__).resolve().parents[3]/"target/debug/modena"))
def base():
 return c.Config(c.Fcc(2.4,[1,1,1]),c.LennardJones(2.5,1.0,1.0),c.VelocityVerlet(0.001),c.Run(steps=4,log_every=1))
class WorkflowTests(unittest.TestCase):
 def test_toml_quotes_and_nested_types(self):
  cfg=base();cfg.structure.element='A"B\\C\nD';cfg.model.types=[{"sigma":1.0,"epsilon":1.0}]
  parsed=tomllib.loads(cfg.to_toml())
  self.assertEqual(parsed["structure"]["element"],cfg.structure.element)
  self.assertEqual(parsed["model"]["types"],[{"sigma":1.0,"epsilon":1.0}])
 def test_typed_surface(self):
  for name in ["NoseHooverChain","Csvr","Ewald","Topology","Constraints","MonteCarlo","Water","RockSalt","Tabulated","Bond","PairException"]:
   self.assertTrue(hasattr(c,name),name)
  cfg=c.Config(c.Fcc(2.4,[1,1,1]),c.LennardJones(2.5,1.0,1.0),c.NoseHooverChain(0.001,1.0),units="reduced",barostat=c.MonteCarlo(0.1,1.0),run=c.Run(steps=3,restart_mode="continue"))
  self.assertIn("OK",check(cfg,BINARY))
 def test_relative_source_paths_and_numeric_rows(self):
  with tempfile.TemporaryDirectory(dir=Path.cwd()) as td:
   source=Path(td)/'a"b.xyz';source.write_text('2\nProperties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1\nA 1 1 1 0 1 0\nB 3 1 1 1 2 0\n')
   cfg=base();cfg.structure=c.Xyz(str(source.relative_to(Path.cwd())),[12.0]*3);cfg.model.types=[{"sigma":1.0,"epsilon":1.0},{"sigma":1.0,"epsilon":0.5}]
   self.assertIn("OK",check(cfg,BINARY));result=run(cfg,BINARY)
   self.assertEqual([r["step"] for r in result.steps_rows],[0,1,2,3,4])
   self.assertTrue(all(isinstance(r["total"],float) and math.isfinite(r["total"]) for r in result.steps_rows))
   self.assertEqual(tomllib.loads(result.toml)["structure"]["file"],str(source.resolve()))
 def test_nonzero_subprocess_propagates(self):
  cfg=base();cfg.integrator.dt=-0.1
  with self.assertRaisesRegex(RuntimeError,"modena run failed"):
   run(cfg,BINARY)

 def test_combined_real_water_and_electrostatic_topology(self):
  water=c.Config(c.Water(8,density=0.3),c.LennardJones(4.0,types=[{"sigma":3.166,"epsilon":0.1554},{"sigma":3.166,"epsilon":0.0}]),c.VelocityVerlet(0.0001),units="real",electrostatics=c.Ewald(0.4,4.0,2.0,method="spme"),thermostat=c.Csvr(100.0,tau=0.05),run=c.Run(steps=3,initial_temperature=100.0,log_every=1))
  self.assertIn("OK",check(water,BINARY));result=run(water,BINARY)
  self.assertEqual(len(result.steps_rows),4);self.assertAlmostEqual(result.steps_rows[0]["temperature"],100.0,places=8)
  with tempfile.TemporaryDirectory() as td:
   src=Path(td)/"dimer.xyz";src.write_text("3\nProperties=species:S:1:pos:R:3:type:I:1:mass:R:1:charge:R:1\nC 2 2 2 0 12 0.1\nH 3.1 2 2 1 1 -0.1\nC 8 5 3 0 12 0\n")
   cfg=c.Config(c.Xyz(str(src),[20.0]*3),c.LennardJones(8.0,types=[{"sigma":3.0,"epsilon":0.2},{"sigma":1.0,"epsilon":0.05}]),c.NoseHooverChain(0.0001,100.0),units="real",electrostatics=c.Ewald(0.4,8.0,2.0),topology=c.Topology(bonds=[c.Bond([0,1],1.1,300.0)],lj_scales=[0.0,0.0,0.5],coulomb_scales=[0.0,0.0,0.8],exceptions=[c.PairException([0,1],0.0,0.0)]),constraints=c.Constraints([c.Constraint([0,1],1.1)]),run=c.Run(steps=3,initial_temperature=100.0,log_every=1))
   self.assertIn("OK",check(cfg,BINARY));result=run(cfg,BINARY)
   self.assertEqual([row["step"] for row in result.steps_rows],[0,1,2,3]);self.assertTrue(all(math.isfinite(row["potential"]) for row in result.steps_rows))
 def test_continuation_state_only_and_intentional_output_paths(self):
  with tempfile.TemporaryDirectory() as td:
   cfg=base();cfg.run.restart_file="resume.chk";cfg.run.restart_every=2
   cfg.run.trajectory_file="traj.xyz";cfg.run.log_file="sentinel.csv"
   sentinel=Path(td)/"sentinel.csv";sentinel.write_text("preserved")
   self.assertIn("OK",check(cfg,BINARY,base_dir=Path(td)));self.assertEqual(sentinel.read_text(),"preserved")
   first=run(cfg,BINARY,base_dir=Path(td));self.assertTrue((Path(td)/"resume.chk").exists());self.assertTrue((Path(td)/"traj.xyz").exists());self.assertEqual(sentinel.read_text(),"preserved")
   cfg.run.restore_file="resume.chk";cfg.run.restart_file="continued.chk"
   second=run(cfg,BINARY,base_dir=Path(td));self.assertEqual([r["step"] for r in second.steps_rows],[4,5,6,7,8]);self.assertAlmostEqual(first.steps_rows[-1]["total"],second.steps_rows[0]["total"],places=9)
   cfg.run.restore_file="continued.chk";cfg.run.restart_mode="state_only";cfg.integrator=c.Langevin(0.001,1.0,1.0)
   self.assertIn("OK",check(cfg,BINARY,base_dir=Path(td)));third=run(cfg,BINARY,base_dir=Path(td));self.assertEqual(third.steps_rows[0]["step"],8)
 def test_nested_dihedral_renderer_and_absolute_table_path(self):
  cfg=base();cfg.topology=c.Topology(dihedrals=[c.Dihedral([0,1,2,3],[c.DihedralTerm(1,0.1,0.0)])])
  parsed=tomllib.loads(cfg.to_toml());self.assertEqual(parsed["topology"]["dihedrals"][0]["terms"][0]["n"],1)
  with tempfile.TemporaryDirectory(dir=Path.cwd()) as td:
   table=Path(td)/'pair"table.dat';table.write_text("0.5 0 0\n1.0 0 0\n2.0 0 0\n3.0 0 0\n4.0 0 0\n")
   cfg=base();cfg.model=c.Tabulated(str(table.relative_to(Path.cwd())))
   self.assertIn("OK",check(cfg,BINARY));result=run(cfg,BINARY)
   self.assertEqual(tomllib.loads(result.toml)["model"]["file"],str(table.resolve()));self.assertTrue(all(r["potential"]==0.0 for r in result.steps_rows))

 def test_toml_del_and_non_ascii_keys_quote_correctly(self):
  cfg=base();cfg.structure.element="label"+chr(127)+"λ"
  parsed=tomllib.loads(cfg.to_toml());self.assertEqual(parsed["structure"]["element"],cfg.structure.element)
