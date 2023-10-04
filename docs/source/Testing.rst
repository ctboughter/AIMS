Testing
=====
As with any good software package, AIMS comes with test data you can use to see if your installation is working properly. All of the examples provided in this documentation utilize the test data in some way, so if you'd like to follow along it is strongly recommended you use the test data.

Thankfully, the pip installation that has become standard as of AIMSv0.9 makes this very easy, and the test data can simply be copied into your current directory (navigated to using the terminal) with the command:

.. code-block:: python

   aims-tests

Importantly, this will not run the tests, but it will copy all of the data into your current directory. From there you should be able to find this test_data directory in either the GUI, the CLI, or the notebook (although the notebook has some code to allow users to run example scripts without copying over the test data).

**Coming Soon**
With the advent of the CLI, it should be very easy and fast to test that every function works. This will happen... eventually.
