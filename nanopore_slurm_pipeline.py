#!/usr/bin/env python3

# To be executed with Python 3.6 or higher
# This queues alignment and variant calling into SLURM handling job's dependencies.

##############################
# Import modules 
##############################

from optparse import OptionParser
import os
import shlex
import subprocess
import re
from configparser import ConfigParser

##############################
# Define functions 
##############################

def create_directory( dir_path ):
  '''
  checks if a directory exists, and if it doesn't, creates it
  '''
  
  if not os.path.exists( dir_path ):
    print( "Creating: %s" %  ( dir_path ))
    os.makedirs( dir_path )
  else:
    print( "DIR already exists: %s\n" % ( dir_path ))

def set_dir_path( top_dir , new_dir ):
  '''
  Takes a path and dirname and joins them intelligently
  '''

  top_dir = os.path.abspath( top_dir )
  return( os.path.join( top_dir , new_dir ))

def setup_workspace( out_dir ):
  '''
  sets up top level nanopore workspace 
  '''

  dir_dict = dict()

  #Setup dict of directory paths 
  dir_dict[ 'output_dir'    ] = out_dir
  dir_dict[ 'alignment_dir' ] = set_dir_path( out_dir , "alignments" )  
  dir_dict[ 'fastq_dir'     ] = set_dir_path( out_dir , "fastq" )   
  dir_dict[ 'variant_dir'   ] = set_dir_path( out_dir , "variant_calling" )
  dir_dict[ 'qc_dir'        ] = set_dir_path( out_dir , "qc" )
  dir_dict[ 'log_dir'       ] = set_dir_path( out_dir , "log" )

  dir_list = [ 'output_dir' , 'alignment_dir', 'fastq_dir', 'variant_dir' , 'qc_dir', 'log_dir' ] 
  for i in dir_list:
    create_directory( dir_dict[ i ] )

  return( dir_dict ) 

def find_fastq( data_dirs ):
  '''
  Finds fastq files in specific directory
  '''

  fastq_list = list()
  for directory in data_dirs:   
   for content in os.listdir( os.path.abspath( directory )):   
    if content.endswith( ".fastq" ) or content.endswith( ".fq.gz" ) or content.endswith( "fastq.gz"):
     fastq_list.append( os.path.join( directory , content ))
  return( fastq_list )

def get_variable( variable, config ):
  '''
  Parses varible from config file
  '''

  parser = ConfigParser()
  with open(config) as stream:
   parser.read_string("[DEFAULT]\n" + stream.read())
  return(parser['DEFAULT'][variable])

def run_command( log_dir , log_pre , script , config, *arg, **kwargs):
  '''
  runs command
  '''
  log = os.path.join( log_dir , log_pre)
  dependency_jobid = kwargs.get('dependency_jobid', None)
  dependency_aligner = kwargs.get('dependency_aligner', None)
                
  if dependency_jobid is not None:
    if dependency_aligner is not None:
      str_command = str("sbatch --error=" +
                log + '.err' +
                " --output=" +
                log + '.out' +
                " --parsable --dependency=afterok:" +
                dependency_jobid + " " +
                script + " "  +
                config + " " +
        dependency_aligner)
    else:
      str_command = str("sbatch --error=" +
                log + '.err' +
                " --output=" +
                log + '.out' +
                " --parsable --dependency=afterok:" +
                dependency_jobid + " " +
                script + " " +
        config)
  else:
    if dependency_aligner is not None:
      str_command = str("sbatch --error=" +
              log + '.err' +
              " --output=" +
              log + '.out' +
              " --parsable " +
              script + " "  +
              config + " " +
          dependency_aligner)
    else:
      str_command = str("sbatch --error=" +
              log + '.err' +
              " --output=" +
              log + '.out' +
              " --parsable " +
              script + " "  +
              config)

  command = subprocess.check_output([str_command], shell=True)
  command = bytes.decode(command)
  command_jobid= re.sub('\n', '', command)
  
  print( "Running " + os.path.basename(log_pre) + " jobid: " + command_jobid )
  print (str_command )
  
  return( command_jobid )

##############################
# Where we do the stuff 
##############################

def main():

  #Setup option parsing 
  parser = OptionParser() 

  parser.add_option( "-c" , "--config" ,  dest="config" ,
                        help="Required: Config file."
                    )
 
  options , args = parser.parse_args()
       
  if not  options.config:
    parser.error( "Please specify a config file.")
  else:
    config = os.path.abspath( options.config ) 
 
  ##Get variables 
  out_dir        = get_variable( 'out_dir', config)
  sample_name    = get_variable( 'sample_name', config )
  scripts_dir    = get_variable( 'scripts_dir', config )
  # cov_source     = get_variable( 'cov_source', config )
  workflow       = get_variable( 'workflow', config )
  data_dirs      = get_variable( 'data_dirs', config)
  data_dirs_list = data_dirs.split(",")
  build = get_variable('build', config)
  build_post = "hg38" if build == "GRCh38" else "hg19"

  ##Define relative paths  
  log_dir         = os.path.join(out_dir, 'log')
  fastq_script    = os.path.join( scripts_dir , 'mergeFastqs.slurm' )
  ngmlr_script    = os.path.join( scripts_dir , 'runNgmlr.slurm' )
  minimap_script  = os.path.join( scripts_dir , 'runMinimap.slurm' )
  sniffles_script = os.path.join( scripts_dir , 'runSniffles.slurm' )
  svim_script = os.path.join( scripts_dir , 'runSvim.slurm' )
  # sniffles_standalone_script = os.path.join( scripts_dir , 'runSnifflesStandalone.slurm' )
  nanosv_script   = os.path.join( scripts_dir , 'runNanoSV.slurm' )
  cov_script      = os.path.join( scripts_dir , 'getCoverage.slurm' )
 
  #step 1 - make workspace
  print( "Making directories" )

  dir_dict = setup_workspace( out_dir ) # Checks if dir exists
  
  # step 2 - make fof from fastq_list

  if not os.path.exists( os.path.join(out_dir, 'fastq', sample_name + '.fastq.fof') ):

      print("Making fof from fastq dir list")
     
      for i in data_dirs_list:
       fastq_list = find_fastq( data_dirs_list )

      fastq_fof = os.path.join(out_dir, 'fastq', sample_name + '.fastq.fof')
      with open(fastq_fof, 'a') as f:
       for item in fastq_list:
        f.write("%s\n" % item)
  else:
    print("Found .fof, skipping .fof creation")

  #step 3 - merge fastq files
  if not os.path.exists( f"{out_dir}/fastq/{sample_name}.merged.fastq" ):
    fastq_jobid = run_command( log_dir = log_dir ,
                             log_pre = 'mergeFastq',
                             script  = fastq_script,
                             config  = config)
  else:
    print("Found merged fastq, skipping merging.")

  fastq_jobid = None

  #step 4 - submit jobs
  for w in workflow.split(','):

    tasks = w.split('-')

    aligner = tasks[0]
    svCallers = tasks[1].split(":") if len(tasks) > 1 else ""

    print("Workflow: " + w)

    if aligner == "minimap":

      if not os.path.exists(f"{out_dir}/alignments/{build}/minimap/{sample_name}.{build_post}.bam"): # We do the alignment, and then the SV calling.
        print(f"Creating minimap2 alignment job for sample {sample_name}")
        if fastq_jobid:
          minimap_jobid = run_command( log_dir = log_dir, 
                     log_pre = 'minimap', 
                     script  = minimap_script, 
                     config  = config,
                     dependency_jobid = fastq_jobid)
        else:
          minimap_jobid = run_command( log_dir = log_dir, 
                     log_pre = 'minimap', 
                     script  = minimap_script, 
                     config  = config)
        if 'sniffles' in svCallers:
          print("Creating Sniffles job")
          sniffles_minimap_jobid = run_command( log_dir    = log_dir, 
                          log_pre    = 'sniffles_minimap',
                          script     = sniffles_script, 
                          config     = config,
                          dependency_jobid = minimap_jobid,
                          dependency_aligner    = 'minimap'
                         )
        if 'svim' in svCallers:
          print("Creating SVIM job")
          svim_minimap_jobid = run_command( log_dir    = log_dir, 
                          log_pre    = 'svim_minimap',
                          script     = svim_script, 
                          config     = config,
                          dependency_jobid = minimap_jobid,
                          dependency_aligner    = 'minimap'
                         )
      else: # Only the SV calling is required
        print("Found minimap alignment, skipping to SV calling")

        if "sniffles" in svCallers:
          print("Creating Sniffles job")
          sniffles_minimap_jobid = run_command( log_dir = log_dir , 
                          log_pre            = 'sniffles_minimap', 
                          script             = sniffles_script, 
                          config             = config, 
                          dependency_aligner = 'minimap'
                         )
        if "svim" in svCallers:
          print("Creating SVIM job")
          svim_minimap_jobid = run_command( log_dir = log_dir , 
                          log_pre            = 'svim_minimap', 
                          script             = svim_script, 
                          config             = config, 
                          dependency_aligner = 'minimap'
                         )
    if aligner == "ngmlr":
      if not os.path.exists(f"{out_dir}/alignments/{build}/ngmlr/{sample_name}.{build_post}.bam"): # We do the alignment, and then the SV calling.
        print(f"Creating NGMLR alignment job for sample {sample_name}")
        if fastq_jobid:
          ngmlr_jobid = run_command( log_dir = log_dir, 
                     log_pre = 'ngmlr', 
                     script  = ngmlr_script, 
                     config  = config,
                     dependency_jobid = fastq_jobid)
        else:
          ngmlr_jobid = run_command( log_dir = log_dir, 
                     log_pre = 'ngmlr', 
                     script  = ngmlr_script, 
                     config  = config)
        if 'sniffles' in svCallers:
          print("Creating Sniffles job")
          sniffles_minimap_jobid = run_command( log_dir    = log_dir, 
                          log_pre    = 'sniffles_ngmlr',
                          script     = sniffles_script, 
                          config     = config,
                          dependency_jobid = ngmlr_jobid,
                          dependency_aligner    = 'ngmlr'
                         )
        if 'svim' in svCallers:
          print("Creating SVIM job")
          svim_minimap_jobid = run_command( log_dir = log_dir, 
                          log_pre    = 'svim_ngmlr',
                          script     = svim_script, 
                          config     = config,
                          dependency_jobid = ngmlr_jobid,
                          dependency_aligner    = 'ngmlr'
                         )
      else: # Only the SV calling is required
        print("Found NGMLR alignment, skipping to SV calling")

        if "sniffles" in svCallers:
          print("Creating Sniffles job")
          sniffles_ngmlr_jobid = run_command( log_dir = log_dir , 
                          log_pre            = 'sniffles_ngmlr', 
                          script             = sniffles_script, 
                          config             = config, 
                          dependency_aligner = 'ngmlr'
                         )
        if "svim" in svCallers:
          print("Creating SVIM job")
          svim_ngmlr_jobid = run_command( log_dir = log_dir , 
                          log_pre            = 'svim_ngmlr', 
                          script             = svim_script, 
                          config             = config, 
                          dependency_aligner = 'ngmlr'
                         )

  #     elif 'nanosv' in variant_caller:
  #       nanosv_ngmlr_jobid = run_command( log_dir               = log_dir , 
  #                                         log_pre               = 'nanosv_ngmlr', 
  #                     script                = nanosv_script, 
  #                         config                = config, 
  #                     dependency_jobid      = ngmlr_jobid, 
  #                     dependency_aligner    = 'ngmlr'
  #                   )    

if __name__ == "__main__":
  import sys
  main()
