
BoreholeNetworksSimulator.unwrap(x::Py) = PythonCall.pyconvert(BoreholeNetworksSimulator.BoreholeOperation, x)

# Account for the different indexing in Python and Julia: Python starts at 0, Julia starts at 1.
BoreholeNetworksSimulator.operate(operator::Py, i, options, X) = operator.operate(i-1, options, X)