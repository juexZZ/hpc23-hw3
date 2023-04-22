# Homework 3 Report

This is the same as Report.md
**Please go to the pdf hw3-report.pdf for the complete homework report!**

## Question 1

(a)

Static scheduling assigns equal size of chunk to the two threads. Therefore, for the first for-loop, the time will be:

    thread 1: 1+2+...+(n-1)/2,

    thread 2: (n-1)/2+1 +  (n-1)/2+2 + ... + n-1

    (let's just assume n is an odd number for simplicity, and this does not affect the generality of our discussion.)

Similarly, for the second for-loop, we have:

    thread 1: n-1 + n-2 + ... + (n-1)/2+1

    thread 2: (n-1)/2 + (n-1)/2-1 + ... + 1

In total, we can see that both of the thread will spend 1 + 2 + ... + n-1 = n(n-1) /2 miliseconds. 

After each for-loop, the threads need to synchronize. For the first for-loop, thread 1 spends less time than thread 2, so it will spend:

[(n-1)/2+1 +  (n-1)/2+2 + ... + n-1] - [1+2+...+(n-1)/2] = (n-1)^2 / 4 

waiting for thread 2.

Similarly, for the second for-loop, the thread 2 will spend the same amount of time (n-1)^2 / 4 waiting for thread 1.

In total, (n-1)^2 / 2 will be spend in waiting (synchronization).

(b): For each for-loop, the execution time will be less imbalanced if we use schedule(static, 1).

For the first for-loop:

    Thread 1: 1 + 3 + 5 + ... + (n-2) = (n-1)^2/4

    Thread 2: 2 + 4 + 6 + ... + (n-1) = (n-1)(n+1) / 4

For the second for-loop, similarly:

    Thread 1: 2 + 4 + 6 + ... + (n-1) = (n-1)(n+1) / 4

    Thread 2: 1 + 3 + 5 + ... + (n-2) = (n-1)^2/4

Combined, the total execution time for both of the thread is the not changed: 1 + 2 + ... + n-1 = n(n-1) /2 miliseconds. But there will be less time spent for waiting to synchronize since it is less imbalanced, so the total running time will be shorter.

(c):

In general, using schedule(dynamic, 1) will improve since it will dynamically assign work to threads and results in better load balancing. But in this specific case, schedule(static, 1) and schedule(dynamic, 1) will be the same. This is because the assignment is not changed given the specific execution time of f(i)s.

(d):

Since f(x) is an independent function, we can use the directive `nowait` for each of the for-loops. This eliminates the implicit synchronization at the end of each for-loop. So at the end each thread will spend exactly n(n-1) /2 miliseconds and no more waiting time.
