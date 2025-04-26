# Maple Transfer Notesheet



Some points I will be addressing in the near future (that you can help start thinking about):

- How can we determine the number of samples we can process at once with our new specs?
  - Generally, optimizing the computing power of our system
- Can we set-up Samantha with a user account that has read & execution permissions so that she can help process data?
- **Spectre:** Needs to be run in conda environment
- **PDMV**: May need to change docker command with RHEL's podman 
- **Sniffles**: Need to finnaly change out of .haplotagged sniffles.
- **ANY DOCKER CONTAINER NEEDS TO BE CLOSED AFTER RUNNING.**
  - John let me know that on Old Maple, there are docker containers still running from >1yr ago. 




## Day 1

**Objective 1:** Set-up and refamiliarize myself with coding / development softwares and workflow.

* Task 1: Fix my terminal
  * Since I won't be ssh'ing on my Mac, Kitty is a strong terminal for nvim editing. I've also configured my nvchad somewhat. Ready to go. 
* Task 2 : Modify `file_transfer.py`
  * Recall the ssh code that transfers files, then modify it to the new PROM ip.
  * Began by 
    * Changing the reference IP to the new prom
    * Added documentation
  * Then added better error checking code
  * And set-up a test environment on PROM to test our scripts
  * SUCCESS!
* Task 3: Clean up `pipeline.py` 
  * Removed all "skipping" options. This was previously used just for bug fixing. Now we don't need it. 
  * Updated the input() functions (that asked step-by-step what files should be used) to one that uses argparse. This has us add all arguments at the command line, and makes it so that we're able to schedule jobs in the future. 

* John Meeting
  * New Maple and Transferred Maple are different things. 
    * We are working with Trasnferred Maple as of now. Soon we will work with New Maple (updated specs)
      *  This is not a problem for today's task
  * On New Maple, we are thinking we don't need a graphical interface. We can use MobaXTerm to launch IGV as a .sh file. 
  * We got old maple data referenceable on New Maple! This means we can open all genomes on IGV!!
    * Now, if we could get firefox on linux, then we could easily upload samples from New Maple too. 
    * We cannot currently write to oldmaple from New Maple. This mainly affect indexing the .vcf files. These will have to be scp'd from oldmaple to a New Maple spot, then reference by IGV for indexing / loading. 



## Day 2:

* Added logging to `file_transfer.py` , tested, validated, and committed ✅
* Added logging to `alignment.py` , tested, validated, and committed ✅
* 

