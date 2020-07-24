RefFinder is a user-friendly web-based comprehensive tool developed  for evaluating and screening reference genes from extensive experimental datasets. It integrates the currently available major computational programs (geNorm, Normfinder, BestKeeper, and the comparative ¦¤Ct method) to compare and rank the tested candidate reference genes. Based on the rankings from each program, It assigns an appropriate weight to an individual gene and calculated the geometric mean of their weights for the overall final ranking. Please cite <a  href="http://www.ncbi.nlm.nih.gov/pubmed/22290409">F Xie, P Xiao, D Chen, L Xu, B Zhang. 2012. miRDeepFinder: a miRNA analysis tool for deep sequencing of plant small RNAs. Plant molecular biology 80 (1), 75-84.</a>

RefFinder is developed by PHP. Users can deploy it to a Php-based server (Apache + PHP) and run it.
If you have a server to host it. Please share your link to other people. 

Thanks!


How to set up RefFinder?

RefFinder was developed in PHP. PHP-supported server should be good for RefFinder, eg. apache+php or XAMPP.

Here is the way how to set up RefFinder to XAMPP on Windows and run it locally:
1. download XAMPP from https://www.apachefriends.org/index.html
  eg. https://www.apachefriends.org/xampp-files/7.4.8/xampp-windows-x64-7.4.8-0-VC15-installer.exe

2. double click xampp-windows-x64-7.4.8-0-VC15-installer.exe to install XAMPP with default settings. 
  XAMPP includes servers (Apache, MySQL, FileZilla FTP server, Mercury Mail Server, and Tomcat), program languages (PHP and Perl), and others.
  
  a. Here only select Apache and PHP
  
  b. Choose a folder to install XAMPP, eg: "D:\XAMPP"
  
  c. Click "Next", "Next", ...
  
  d. Finally, click "Finish" and launch XAMPP Control panel (XAMPP Control panel can be also opened by clicking D:\XAMPP\XAMPP-control.exe)
  
  e. Start Apache (Click Apache's "Start" and you will see the ports like "80, 443")

3. Download RefFinder's source code from https://github.com/fulxie/RefFinder/archive/master.zip
4. Unzip RefFinder-master.zip and copy RefFinder-master folder to "htdocs" of XAMPP installation folder (eg: D:\XAMPP\htdocs\RefFinder-master). 
5. Open Internet browser and type https://localhost/RefFinder-master/index.php in address bar, then you can use RefFinder now.
6. In future, you need to use XAMPP control panel (D:\XAMPP\XAMPP-control.exe) to start Apache first and then type the url (https://localhost/RefFinder-master/index.php) in the browser to use RefFinder


