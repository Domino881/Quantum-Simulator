# Quantum-Simulator
Hobby project - creating a simple quantum circuit simulator in C++.
Inspired by https://github.com/qiskit/

Example output:
```
|0> --[H]---*--------------------------------------------------------------------------I----x---[H]---*---------*--------------*---[M]----------------
            |                                                                          I    |         |         |              |
|0> --[H]---|----*----*----------------------------------------------------------------I----|----x---[P]--[H]---|----*---------|----*---[M]-----------
            |    |    |                                                                I    |    |              |    |         |    |
|0> --[H]---|----|----|----*----*----*----*--------------------------------------------I----|----x-------------[P]--[P]--[H]---|----|----*---[M]------
            |    |    |    |    |    |    |                                            I    |                                  |    |    |
|0> --[H]---|----|----|----|----|----|----|----*----*----*----*----*----*----*----*----I----x---------------------------------[P]--[P]--[P]--[H]--[M]-
            |    |    |    |    |    |    |    |    |    |    |    |    |    |    |    I
|0> --[X]--[P]--[P]--[P]--[P]--[P]--[P]--[P]--[P]--[P]--[P]--[P]--[P]--[P]--[P]--[P]---I--------------------------------------------------------------

00000: 5
00001: 5
00010: 6
00011: 14
00100:# 37
00101:################################## 690
00110:######### 187
00111:# 28
01000: 2
01001: 10
01010: 1
01011: 4
01100: 5
01101: 3
01110: 2
01111: 1
```
