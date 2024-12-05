# Chrono Concrete

&#x20;[![BSD License](http://www.projectchrono.org/assets/logos/chrono-bsd.svg)](LICENSE)

This is a GitLab repository for the development and implementation of concrete modeling codes in the [Project Chrono](https://www.projectchrono.org) software. The repository and wiki are in development and may be incomplete or include mistakes.

## Installation and Setup

Select below instructions for installation and setup on Northwestern's Quest HPC system or local desktop. Note that Quest is the preferred platform and desktop instructions are incomplete at this time.

### ☁ Quest Installation and Setup

Follow the below instructions to get set up with Project Chrono on Quest. These instructions assume that you have a Northwestern **NetID** and access to the group project allocation **p31861**. If you do not have either of these, please contact a Chrono Concrete admin.

Note 1: You may also have access to the Northwestern group Buy-in node **b1187**. This node can be used the same as with the generic project allocation; just swap p31861 for b1187. Even if using the buy-in node, you should still store files within the project directory though. Additionally, remember that if running on the buy-in node, you will be limited to requesting the number of processors that are physically located on that node.

Note 2: When following these instructions, if you encounter errors, please try manually typing in the commands rather than copying them from your browser window.

<details>

<summary>Step 1: Install and Configure SSH Client </summary>

Install an SSH Client

* Any SSH client should work, but we recommend PuTTY
* Download and install PuTTY from [https://www.putty.org/](https://www.putty.org/)
* After installation open PuTTY. Within the PuTTY Configuration window:
  * Enter '**quest.northwestern.edu**' for Host Name
  * Enter '**22**' for Port&#x20;
  * Enter '**Quest**' for Saved Sessions
  * Click **Save**
  * Click **Quest** which should not be added to the Session list
  * Click **Open**
* A new SSH window will open. In this window login with your Northwestern NetID and password

</details>

<details>

<summary>Step 2: Install and Configure SFTP Client</summary>

Install an SFTP Client

* Any SFTP client should work, but we recommend FileZilla
* Download and install FileZilla Client from [https://filezilla-project.org/](https://filezilla-project.org/)
* After installation open FileZilla. Within the FileZilla window:
  * Click **File** and then **Site Manager...**
  * In the opened window click **New Site** and enter 'Quest' for the name
  * Enter '**quest.northwestern.edu**' for Host&#x20;
  * Enter '**22**' for Port&#x20;
  * Enter your NetID for **User** and password for **Password**
  * Click **New Bookmark** and enter 'Projects' for the name
  * Choose any Local directory that you want
  * Enter '**/projects/p31861**' for Remote directory
  * Click **OK**
  * Click **File** and then **Site Manager...**
  * Click **Connect**
* The remote site on the right side of your window should automatically connect to the Quest Project Chrono Project and you should see a folder called **Singularity Container**

</details>

<details>

<summary>Step 3: Create Your User Directory</summary>

Create a directory for all of your developments and testing. No files/folders should be created or changed at the top-most '/projects/p31861' directory.

* Within FileZilla, enter the **Users** folder&#x20;
* Right click in the '/projects/p31861/Users' folder and select **Create Directory and Enter It**
* Name the folder with your name in the following format **LastnameFirstname**

</details>

<details>

<summary>Step 4: Copy Singularity File to User Directory</summary>

Copy the SIF file into your User Directory

* In the SSH window run the following command, being sure to replace **LastnameFirstname** with your correct directory name

<pre><code><strong>cp /projects/p31861/SingularityContainer/project-chrono-dependencies-with-intel-mkl.sif /projects/p31861/Users/LastnameFirstname 
</strong></code></pre>

</details>

<details>

<summary>Step 5: Clone Chrono-Concrete</summary>

Clone the Project Chrono GitHub into your User Directory

* In the SSH window cd into your User Directory with the following command, being sure to replace **LastnameFirstname** with your correct directory name&#x20;

<pre><code><strong>cd /projects/p31861/Users/LastnameFirstname
</strong></code></pre>

* Clone the GitHub project here with the following command

```
git clone https://github.com/Concrete-Chrono-Development/chrono-concrete.git
```

* Pull updates to GitHub project - Make sure to manually type in these commands or they may not work when copy-pasted.

```
cd chrono-concrete
git pull https://github.com/Concrete-Chrono-Development/chrono-concrete.git
git submodule init​
git submodule update
```

</details>

<details>

<summary>Step 6: Build Chrono-Concrete</summary>

Copy example make script, edit, and build Project Chrono

* Copy example make script to User Directory, being sure to replace **LastnameFirstname** with your correct directory name&#x20;

<pre><code><strong>cp /projects/p31861/ExampleScripts/submit_chrono_make.sh /projects/p31861/Users/LastnameFirstname 
</strong></code></pre>

* Navigate to the newly copied 'submit\_chrono\_make.sh' file in FileZilla and double-click on it to edit. Change all instances of **LastnameFirstname** in the file to your appropriate directory and save/upload editted file back to Quest
* In your SSH client navigate to your User Directory and run the following command to submit job

```
sbatch submit_chrono_make.sh
```

You can check the status of your job with the command, being sure to replace **NetID** with your NetID:

```
squeue -u NetID
```

Once the job has completed, proceed to Step 7.

</details>

<details>

<summary>Step 7: Verify Installation</summary>

Verify proper installation of Chrono-Concrete by running a test job with MPI

* Copy example make script to User Directory, being sure to replace **LastnameFirstname** with your correct directory name&#x20;

```
cp /projects/p31861/ExampleScripts/example_submit_mpi.sh /projects/p31861/Users/LastnameFirstname 
```

* Navigate to the newly copied 'submit\_chrono\_make.sh' file in FileZilla and double-click on it to edit. Change all instances of **LastnameFirstname** in the file to your appropriate directory and save/upload editted file back to Quest
* Make an output directory, being sure to replace **LastnameFirstname** with your correct directory name&#x20;

```
mkdir /projects/p31861/Users/LastnameFirstname/outdir
```

* In your SSH client navigate to your User Directory and run the following command to submit job

```
sbatch example_submit_mpi.sh
```

You can check the status of your job with the command, being sure to replace **NetID** with your NetID:

```
squeue -u NetID
```

Once the job has completed, open the outlog file in your User Directory and confirm that the simulation ran. Then navigate to the output directory (./outdir/TestJob) and confirm that several .csv files were created.&#x20;

</details>

<details>

<summary>Step 8: Development</summary>

Code within the chrono-concrete directory can be developed as needed and be pushed/pulled to the GitHub. Please read online about how git works so that you ensure you are properly developing with everyone else.&#x20;

You can modify/copy the example .sh scripts and outdir in your User Directory to help your developments.

Please **do not** edit any files outside of your User Directory.

</details>


### 💻 Desktop Installation and Setup

Follow the below instructions to get set up with Project Chrono. These instructions were written with Windows users in mind. It is also possible to install on Mac or Linux; Mac instructions are located in the notes section of each Step.

Please note the versions of each piece of software. Newer or alternate versions may work, but have not been tested and verified for compatibility.

<details>

<summary>Step 1: Install C++ Compiler</summary>

Install Microsoft Visual Studio 2022. The [community edition of the latest Visual Studio](https://visualstudio.microsoft.com/downloads/) is available for free.

* During installation make sure to check and include _"Desktop Development with C++"_
* After installation open Visual Studio and sign-in if necessary

```
- For Mac: Use Xcode Package. Download via App Store for free - it contains the clang++ compiler.
- Notes: Other compilers were also tested (e.g. Intel C++, PGI) but they are not officially supported and maintained. While it is likely possible to build Chrono with other toolchains, this might require changes to the CMake scripts.
- Notes: Any version of Visual Studio 2019 or newer should work. Visual Studio 2017 has problems with the heavy use of inlining in recent version of Eigen. For the latest version of Chrono (specifically due to the reimplemented ANCF elements), this can result in very long compilation times or even hang ups. We recommend using VS 2019 or newer.
```

</details>

<details>

- C++ API reference
  - [Main development branch](http://api.projectchrono.org/)
  - [Release 9.0.0](http://api.projectchrono.org/9.0.0/)
  - [Release 8.0.0](http://api.projectchrono.org/8.0.0/)
- Python interface
  - [PyChrono](https://api.projectchrono.org/pychrono_introduction.html)
- Reference manuals
  - [Core module](https://api.projectchrono.org/manual_root.html)
  - [Chrono::Vehicle module](https://api.projectchrono.org/manual_vehicle.html)

Download the Eigen version 3.4.0 zipped  source code. This [code is available for free](https://gitlab.com/libeigen/eigen/-/releases/3.4.0).&#x20;

* Unzip the downloaded file and store the contents in an easy-to-find location, suggested location is: C:\workspace\libraries\eigen-3.4.0

<pre><code><strong>- For Mac: Install it via homebrew: brew install eigen. Homebrew installs into /opt/homebrew since MacOS 12 Monterey and the new Apple Silicon (arm46, M1, M2...) hardware. If Eigen is not found automatically, you can search its folder with:
</strong>find /opt/homebrew -name Eigen
<strong>- Notes: Any version of Eigen 3.3.0 or newer should work. 
</strong></code></pre>

</details>

<details>

<summary>Step 3: Install Irrlicht Library</summary>

Download Irrlicht SDK version 1.8.5. This [code is available for free](https://irrlicht.sourceforge.io/?page\_id=10).&#x20;

* Unzip the downloaded file and store the contents in an easy-to-find location, suggested location is: C:\workspace\libraries\irrlicht-1.8.5

<pre><code><strong>- For Mac: The best way to install irrlicht on the Mac is: brew install irrlicht (release v.1.8.5). On MacOS 12 (Monterey) you have to set IRRLICHT_ROOT to /opt/homebrew.
</strong><strong>- Notes: Any version of Irrlicht SDK version 1.8.2 or newer should work. 
</strong></code></pre>

</details>

<details>

<summary>Step 4: Install CMake</summary>

Install CMake version 3.25.0. An installer for the [software is available for free](https://cmake.org/download/).

* During installation be sure to check for the CMake executable to be included in your Path environmental variable

<pre><code><strong>- For Mac: The CMake.app bundle also contains command line tools, you must set appropriate links to use it from the terminal. It is better to install a pure command line version via homebrew (https://brew.sh). After installing the home brew package manager type: brew install cmake in the terminal.
</strong><strong>- Notes: Any version of CMake version 1.8.2 or newer should work. 
</strong></code></pre>

</details>

<details>

<summary>Step 5: Install Git Client - In Progress</summary>

*

</details>

<details>

<summary>Step 6: Clone Project Chrono - <em>In Progress</em></summary>

*

</details>

<details>

<summary>Step 7: Run CMake - <em>In Progress</em></summary>

*

</details>

<details>

<summary>Step 8: Compile Project - <em>In Progress</em></summary>

*

</details>

Source for many of these instructions is Project Chrono documentation. More details [here.](https://api.projectchrono.org/tutorial\_install\_chrono.html)

## Questions

All questions and feature requests (related to this fork of Project Chrono) should be directed to the Discussions at the top of the Repository. Questions related to Project Chrono in general should be directed to their own forum here: https://groups.google.com/g/projectchrono

> #### [Discussions](https://github.com/Concrete-Chrono-Development/chrono-concrete/discussions)

## Issues and Bugs

All found issues and bugs should be directed to the Issues at the top of the Repository.

> #### [Issues](https://github.com/Concrete-Chrono-Development/chrono-concrete/issues)

## Wiki, Details, and Examples

All examples, details, and examples are available on the GitHub Wiki which is located at the Wiki tabe at the top of this page, or you can access it at the link below:

> #### [Chrono Concrete Wiki](https://github.com/Concrete-Chrono-Development/chrono-concrete/wiki)

## Error Fixes

A list of common errors and their fixes are available on the Wiki at the below page. If you encounter an error and find a solution to it, please include it on this page.

> #### [Wiki Error Fixes](https://github.com/Concrete-Chrono-Development/chrono-concrete/wiki)

## License

Project Chrono, and by extension the Chrono Concrete application, is licensed free and open source under the [BSD 3-Clause “New” or “Revised” License](https://choosealicense.com/licenses/bsd-3-clause/). One reason for Project Chrono's popularity is that its users are granted the freedom to modify and redistribute the software and have a right of continued free use, within the terms of the BSD license.

## Acknowledgements

This research was supported in part through the computational resources and staff contributions provided for the Quest high performance computing facility at Northwestern University which is jointly supported by the Office of the Provost, the Office for Research, and Northwestern University Information Technology.

## Citing this Work

A proper citation for theses developments is still to come. In the meantime, please reference this GitHub in all publications.
