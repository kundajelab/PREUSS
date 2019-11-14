# -*- coding: utf-8 -*-

"""
Python Object serialize/deserialize Solution
================

It aims to help decode / encode "JSON-incompatible" python object so that it can be passed as `string`.

Highlights of the solution:
- It is achieved by overdrive the `JSONEncoder.default()`.
- For serialization/deserialization, it is done by `pickle`.

Ref - https://code.tutsplus.com/tutorials/serialization-and-deserialization-of-python-objects-part-1--cms-26183

"""

from typing import Dict, Any

import json
from json import JSONEncoder
import pickle

import base64
import numpy as np

# ----------------------------------
# region Python Object

PYTHON_OBJECT_JSON_KEY = '_python_object'


class PythonObjectEncoder(JSONEncoder):
    r"""
    Custom JSON Encoder - Python Object Encoder

    It utilizes `pickle` to help do the work.

    NOTE:
        - In Python 3.x, `json.dumps()` returns `bytes` instead of `string`, which `JSONEncoder` cannot handle.
        - To work around, "decode" the dumps by `latin1` and "encode" back when "loading" it.
            - This works because arbitrary binary strings are valid latin1 which can always be decoded to Unicode
                and then encoded back to the original string again.
        - Ref - https://stackoverflow.com/questions/18478287/making-object-json-serializable-with-regular-encoder

    """

    def default(self, python_obj: Any) -> Dict[str, Any]:
        return {
            # Use `HIGHEST_PROTOCOL` to reduce the encoding size.
            PYTHON_OBJECT_JSON_KEY: pickle.dumps(python_obj, protocol=pickle.HIGHEST_PROTOCOL).decode('latin1')
        }


def as_python_object(serialized_dict: Dict[str, Any]) -> Any:
    r"""
    "Object Hook" function for `json.loads()`.

    It help load the serialized Python object by custom `PythonObjectEncoder`.

    Parameters
    ----------
    serialized_dict: Dict[str, Any]
        Passed serialized data dict.

    Returns
    -------
    python_object: Any
        The deserialized Python object.
    """

    if PYTHON_OBJECT_JSON_KEY in serialized_dict.keys():
        return pickle.loads(serialized_dict[PYTHON_OBJECT_JSON_KEY].encode('latin1'))

    # Default
    return {}

# endregion


# ----------------------------------
# region 2-D Data Array

def encode_base64(ndarray):
    r"""
    Encode a "numpy array" data via Base64.

    Parameters
    ----------
    ndarray: np.array
        The 2-D array, in "numpy array" format

    Returns
    -------
    json_dump: List
        The output JSON dump, a list of elements - "data type", "base64-encoded data", "shape".

    """
    return json.dumps([str(ndarray.dtype), base64.b64encode(ndarray), ndarray.shape])


def decode_base64(json_dump):
    r"""
    Decode a JSON dump data via Base64.

    Parameters
    ----------
    json_dump: List
        JSON dump data, a list of elements - "data type", "base64-encoded data", "shape".

    Returns
    -------
    numpy_arr: numpy.array
        The 2-D array, in "numpy array" format

    """
    loaded = json.loads(json_dump)
    dtype = np.dtype(loaded[0])
    arr = np.frombuffer(base64.decodestring(loaded[1]), dtype)

    # Check if need to do "re-shape"
    if len(loaded) > 2:
        return arr.reshape(loaded[2])
    return arr

# endregion
