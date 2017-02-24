import unittest,os
from glob import glob

class TestCCode(unittest.TestCase):

  def setUp(self):
    self.toTest = ["./lungs_analyze",\
    "/usr/local/share/projects/COPD/Data/EPOC_patients/ENFISEMA 01/1.3.12.2.1107.5.1.4.60077.30000007111412001959300000060/"]
    
  def testBuild(self):
    if not os.path.exists('lungs_analyze'):
	self.fail("No bin to test!!")
    print "Testing ",self.toTest
  
    retcode = os.spawnlp(os.P_WAIT, self.toTest[0], self.toTest[0], self.toTest[1])
    self.failUnless(retcode==0,"Program exited with the retcode"+str(retcode))
    
    files = glob("./*.dcm")
    print files,len(files)

def main():
  unittest.main()

if __name__ == '__main__':
  main()
