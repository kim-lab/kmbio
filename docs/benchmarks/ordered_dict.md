# Use `OrderedDict`

## Entire test suite (without mmCIF vs. mmCIF ref.)

### Using `_child_list` and `_child_dict`

-   Time: 14.5 minutes
-   Memory: 4.5 GB

```text
Command being timed: "py.test"
User time (seconds): 790.59
System time (seconds): 8.48
Percent of CPU this job got: 91%
Elapsed (wall clock) time (h:mm:ss or m:ss): 14:28.96
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 4470944
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 151
Minor (reclaiming a frame) page faults: 2204043
Voluntary context switches: 56095
Involuntary context switches: 32386
Swaps: 0
File system inputs: 109672
File system outputs: 1331104
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
```

### Using `OrderedDict` (`_children`)

-   Time: 14 minutes
-   Memory: 3.7 GB

```text
Command being timed: "py.test"
User time (seconds): 766.70
System time (seconds): 6.15
Percent of CPU this job got: 91%
Elapsed (wall clock) time (h:mm:ss or m:ss): 14:02.37
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 3696292
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 31
Minor (reclaiming a frame) page faults: 1851529
Voluntary context switches: 62037
Involuntary context switches: 25256
Swaps: 0
File system inputs: 60232
File system outputs: 1331096
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
```

## Test five large mmCIF files

### Using `_child_list` and `_child_dict`

-   Time: 7 minutes 19 seconds
-   Memory: 1.56 GB

```text
Command being timed: "py.test tests/test_bioassembly.py::test_mmcif_vs_mmcif_ref"
User time (seconds): 427.67
System time (seconds): 4.05
Percent of CPU this job got: 98%
Elapsed (wall clock) time (h:mm:ss or m:ss): 7:19.94
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 1558948
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 95
Minor (reclaiming a frame) page faults: 1100988
Voluntary context switches: 14205
Involuntary context switches: 37223
Swaps: 0
File system inputs: 6136
File system outputs: 383560
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
```

### Using `OrderedDict`

-   Time: 7 minutes 11 seconds
-   Memory: 1.30 GB

```text
Command being timed: "py.test tests/test_bioassembly.py::test_mmcif_vs_mmcif_ref"
User time (seconds): 420.45
System time (seconds): 3.04
Percent of CPU this job got: 98%
Elapsed (wall clock) time (h:mm:ss or m:ss): 7:11.05
Average shared text size (kbytes): 0
Average unshared data size (kbytes): 0
Average stack size (kbytes): 0
Average total size (kbytes): 0
Maximum resident set size (kbytes): 1293356
Average resident set size (kbytes): 0
Major (requiring I/O) page faults: 15
Minor (reclaiming a frame) page faults: 920721
Voluntary context switches: 10203
Involuntary context switches: 167275
Swaps: 0
File system inputs: 23504
File system outputs: 383336
Socket messages sent: 0
Socket messages received: 0
Signals delivered: 0
Page size (bytes): 4096
Exit status: 0
```
