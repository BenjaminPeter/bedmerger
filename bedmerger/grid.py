import sys, socket

"""
runs the merging script on the sun grid
"""



def merge( params ):
    """
    """

    pass


    python_exe = sys.executable

    setup_script="""
#!/bin/bash
#$ -N test
# Giving the name of the output log file
#$ -o test.log
# Combining output/error messages into one file
#$ -j y
# One needs to tell the queue system to use the current directory as the working directory
# Or else the script may fail as it will execute in your top level home directory /home/username
#$ -cwd
# Now comes the commands to be executed

python 

/share/apps/matlab/bin/matlab -nodisplay -nodesktop -nojvm -r matlab-test
# Note after -r is not the name of the m-file but the name of the routine
exit 0
"""
