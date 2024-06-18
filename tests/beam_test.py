# test_beam_module.py

from beam_module import Beam

# TODO: test that beam can be properly instantiated (all attributes are what they should be)
def test_class_attr():
    beam = Beam()
    assert beam.t_0 == 1 # or whatever we set it to in a provided txt file