import importlib.util, tempfile, unittest
from pathlib import Path

spec=importlib.util.spec_from_file_location("integration",Path(__file__).parents[1]/"bin"/"integrate_regulatory.py")
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)

class IntegrationUnitTests(unittest.TestCase):
    def test_half_open_overlap(self):
        self.assertEqual(m.overlap(0,10,10,20),0)
        self.assertEqual(m.overlap(0,11,10,20),1)
    def test_version_harmonization(self):
        self.assertEqual(m.base_gene("ENSG00000123456.7"),"ENSG00000123456")
        self.assertEqual(m.base_gene("SYMBOL"),"SYMBOL")
    def test_states(self):
        self.assertEqual(m.state(1,.01,True,.05),"UP")
        self.assertEqual(m.state(-1,.01,True,.05,True),"HYPOMETHYLATED")
        self.assertEqual(m.state(1,.2,True,.05),"NO_SIGNIFICANT_CHANGE")
        self.assertEqual(m.state(1,.01,False,.05),"NOT_TESTED")
    def test_interval_rejection(self):
        with self.assertRaises(ValueError): m.validate_interval_rows([{"region_id":"x","chromosome":"chr1","start":"10","end":"10"}],"test")
        with self.assertRaises(ValueError): m.validate_interval_rows([{"region_id":"x","chromosome":"chrM","start":"0","end":"1"}],"test")

if __name__=="__main__": unittest.main()
