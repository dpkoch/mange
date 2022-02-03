import sys

import numpy as np
import pytest

from mange_python import SE2, SE3, SO2, SO3

GROUP_TYPES = [SE2, SE3, SO2, SO3]


@pytest.mark.parametrize("G", GROUP_TYPES)
def test_constructor(G):
    # Simply make sure this doesn't raise an exception or segfault
    G()


@pytest.mark.parametrize("G", GROUP_TYPES)
def test_exp(G):
    assert isinstance(G.Exp(np.random.random(size=(G.DOF,))), G)


@pytest.mark.parametrize("G", GROUP_TYPES)
def test_log(G):
    if G == SO2:
        assert isinstance(G.Random().Log(), float)
    else:
        assert G.Random().Log().shape == (G.DOF,)


@pytest.mark.parametrize("G", GROUP_TYPES)
def test_hat_vee_roundtrip(G):
    vector = np.random.random(size=(G.DOF,))
    assert np.all(G.vee(G.hat(vector)) == vector)


@pytest.mark.parametrize("G", GROUP_TYPES)
def test_self_multiply(G):
    assert isinstance(G.Random() * G.Random(), G)


@pytest.mark.parametrize("G", GROUP_TYPES)
def test_vector_multiply(G):
    vector = np.random.random(size=(G.DIM,))

    assert (G.Random() * vector).shape == (G.DIM,)
    assert np.all(G.Identity() * vector == vector)


if __name__ == "__main__":
    sys.exit(pytest.main([__file__]))
