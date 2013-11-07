import subprocess

nSessions = 5

for i in xrange(nSessions):
  subprocess.Popen(["./spatiocyte-core"])
