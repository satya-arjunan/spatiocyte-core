import subprocess

nSessions = 3

for i in xrange(nSessions):
  subprocess.Popen(["./spatiocyte-core"])
