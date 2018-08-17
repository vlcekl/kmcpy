import pytest
sys.path.append("../../../kmcsim/")
import kmcsim

def test_make_fcc(box):
    assert box.shape[0] == 3, "make_fcc(box): Wrong box shape!"
