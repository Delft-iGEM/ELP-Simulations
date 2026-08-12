import os

from calvados.cfg import Config, Components
from calvados.sim import Sim
from pathlib import Path

from tools.elibpy import build_sequence_with_features

# Parameters
# OpenMM runs in nm where as PDB is in Angtrums so divide it by 10
z_wall = 230 / 10

def build_sim(sim: Sim):
     pass

# Job settings for Delft Blue
partition = "gpu-a100-small"
runtime = "3:30:00"
cpu_per_task = "2"

# or
# partition = "gpu-a100"
# runtime = "24:30:00"
# cpu_per_task = "18"

sim_name = Path(__file__).parent.name

L = 50
N_save = 7000
N_frames = 300

'''
sequences: dict[str, str] = {
     "VPGIG70": build_sequence_with_features([
          ("K70", "K")
     ])["seq"]
}


sequences: dict[str, str] = {
     "VPGIG70": build_sequence_with_features([
          ("K18", "K"),
          ("I34", "I"),
          ("K18", "K")
     ])["seq"]
}
'''
sequences: dict[str, str] = {
     "VPGIG70RGD": build_sequence_with_features([
          ("I18", "I"),
          ("K17", "K"),
          ("I17", "I"),
          ("K18", "K"),
          ("RGD","RGD")
     ])["seq"]
}

sequences["VPGIG70RGD"] = f"Z{sequences['VPGIG70RGD'][1:]}"

if __name__ == "__main__":
     path = Path(__file__).parent.resolve()
     cwd = Path(os.getcwd())

     runtime_dir = path / "runtime"
     runtime_dir.mkdir(exist_ok=True)

     fasta_file = runtime_dir / "molecules.fasta"

     # Config preparation
     config = Config(
          sysname = sim_name,
          box = [L, L, L],
          temp = 293.15,
          ionic = 0.13,
          pH = 7.5,
          ext_force = False,
          ext_force_expr = f'step({z_wall}-z)*0.5*({z_wall}-z)^2',
          topol = 'center',   
          wfreq = N_save,
          steps = N_frames*N_save,
          runtime = 0,
          platform = 'OpenCL',
          restart = 'checkpoint',
          frestart = 'restart.chk',
          verbose = True
     )

     components = Components(
          # Defaults
          molecule_type = 'protein',
          nmol = 1, # number of molecules
          restraint = False,
          charge_termini = 'both',
          fresidues = str(cwd / "residues_CALVADOS2.csv"),
          ffasta = str(fasta_file)
     )

     for k in sequences.keys():
          components.add(name=k)

     # Prepare runtime - You probably don't need to edit below this

     # Write config and components
     config.write(str(runtime_dir), "config.yaml")
     components.write(str(runtime_dir), "components.yaml")

     # Create FASTA
     fasta_content = "\n".join([f">{k}\n{v}" for k, v in sequences.items()])
     
     fasta_file.write_text(fasta_content)
     
     # Move job.sh and run.py
     run_file = (cwd / "template/run.py").read_text()
     (runtime_dir / "run.py").write_text(run_file)

     job_file = (cwd / "template/job.sh").read_text()
     
     job_file = job_file \
          .replace("{sim_name}", sim_name) \
          .replace("{partition}", partition) \
          .replace("{runtime}", runtime) \
          .replace("{cpu_per_task}", cpu_per_task)

     (runtime_dir / "job.sh").write_text(job_file)