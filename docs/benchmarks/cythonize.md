# Cythonize `mmcif2dict`

## Entire test suite

```text
Command being timed: "py.test"
User time (seconds): 241.72
System time (seconds): 6.41
Percent of CPU this job got: 81%
Elapsed (wall clock) time (h:mm:ss or m:ss): 5:03.98
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 3882396
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 0
Minor (reclaiming a frame) page faults: 2153011
Voluntary context switches: 31231
Involuntary context switches: 23344
Swaps: 0
File system inputs: 40
File system outputs: 1330768
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
```

## Test 5 large mmCIF files

```text
Command being timed: "py.test tests/test_bioassembly.py::test_mmcif_vs_mmcif_ref"
User time (seconds): 173.65
System time (seconds): 3.25
Percent of CPU this job got: 95%
Elapsed (wall clock) time (h:mm:ss or m:ss): 3:05.38
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 1286308
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 0
Minor (reclaiming a frame) page faults: 946342
Voluntary context switches: 10670
Involuntary context switches: 15304
Swaps: 0
File system inputs: 104
File system outputs: 383384
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
```
