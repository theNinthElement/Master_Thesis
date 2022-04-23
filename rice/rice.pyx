cdef extern from "rice.h":
    int RiceEncodeFile(FILE *inFile, FILE *outFile, const unsigned char k);
    int RiceDecodeFile(FILE *inFile, FILE *outFile, const unsigned char k);

def py_RiceEncodeFile(inFile: bytes, outFile: bytes, k: bytes) -> None:
    RiceEncodeFile(inFile, outFile, k)

def py_RiceDecodeFile(inFile: bytes, outFile: bytes, k: bytes) -> None:
    RiceDecodeFile(inFile, outFile, k)
