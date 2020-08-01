<?php


 if(isset($_GET['sessionID'])) {
        $sessionID=$_GET['sessionID'];
    }

if(isset($_GET['type']) ? $_GET['type'] : null) {

         $type = $_GET['type'];
        if($type == 'ssea' || $type == 'msea')
        {

            $genefile = "./Data/Pipeline/Resources/shinyapp2_temp/$sessionID".".SSEA2PHARM_genes.txt";
            $resultfile="./Data/Pipeline/Results/shinyapp2/$sessionID".".SSEA2PHARM_app2result.txt"; 

         }
         else
         {
          $genefile = "./Data/Pipeline/Resources/shinyapp2_temp/$sessionID".".KDA2PHARM_genes.txt";
          $resultfile="./Data/Pipeline/Results/shinyapp2/$sessionID".".KDA2PHARM_app2result.txt"; 
         }   

         if($type == 'pharm')
         {
            $resultfile="./Data/Pipeline/Results/shinyapp2/$sessionID"."_app2result.txt"; 
         }
    }
    


function scientificNotation($val){
    $exp = floor(log($val, 10));
    if($val == 0)
    {
      return 0;
    }
    else
    {
       return sprintf('%.3fE%+03d', $val/pow(10,$exp), $exp);
    }
   
}




if (!file_exists($resultfile)) 
{

    shell_exec("sh run_app2.sh $sessionID | tee ./Data/Pipeline/Results/shinyapp2/out.txt");
    $outfile = "./Data/Pipeline/Results/shinyapp2/out.txt";
    sleep(1);
    if (file_exists($outfile)) {
        unlink($outfile);
    } 
     chmod($resultfile, 0777);
     

    
}


 ?>


  <table class="table table-bordered review" style="text-align: center"; id="shinyapp2resultstable">
<thead>
     <tr >
        <th colspan="2">
          Download Output Files
        </th>
      </tr>
      <tr>
        <td>
          Network Based Drug Repositioning Analysis File
        </td>
        <td>
          <a href=<?php print($resultfile); ?> download> Download</a>
        </td>
      </tr>
      
      <tr>

</thead>
</table>



  <table id="shinyapp2_results" class="table table-striped table-bordered" cellspacing="0" width="100%">
    <thead>
      <tr>
        <th>
            Drug Name
        </th>
        <th>
            Z-Score
        </th>
        <th>
            Jaccard Score
        </th>

      </tr>
    </thead>
    <?php 
        $result = file($resultfile);
        foreach($result as $line)
        {
          $line_array = explode("\t", $line);
          $drug = trim($line_array[0]);
          $zscore = number_format($line_array[1], 3, ".", "");
          $jaccard = scientificNotation(trim((float)$line_array[2]));
          echo "<tr><td>$drug</td><td>$zscore</td><td>$jaccard</td></tr>";
        }

    ?>
  </table>

<?php

$results_sent="./Data/Pipeline/Results/shinyapp2_email/$sessionID"."sent_results";
$email = "./Data/Pipeline/Results/shinyapp2_email/$sessionID"."email";

if ((!(file_exists($results_sent)))) 
{
    if(file_exists($email)) 
    {
        require('./PHPMailer-master/class.phpmailer.php');

        $resultfile = "./Data/Pipeline/Results/shinyapp2/$sessionID"."_app2result.txt";
        
        $emailid="./Data/Pipeline/Results/shinyapp2_email/$sessionID"."email";
        $mail = new PHPMailer();

        $mail->Body = 'Congratulations! You have successfully executed our pipeline. Please download your results.';
        $mail->Body .= "\n";
        $mail->Body .= 'Your results are available at: http://mergeomics.research.idre.ucla.edu/runpharmomics.php?sessionID=';
        $mail->Body .= "$sessionID";
        $mail->Body .= "\n";
        $mail->Body .= 'Note: Your results will be deleted from the server after 24 hours';

        //$mail->IsSMTP(); // telling the class to use SMTP

        $mail->SMTPAuth   = true;                  // enable SMTP authentication
        $mail->SMTPSecure = "tls";                 // sets the prefix to the servier
        $mail->Host       = "smtp.gmail.com";      // sets GMAIL as the SMTP server
        $mail->Port       = 587;                   // set the SMTP port for the GMAIL server
        $mail->Username   = "dougvarneson@gmail.com";  // GMAIL username
        $mail->Password   = "friday11";            // GMAIL password

        $mail->SetFrom('darneson@g.ucla.edu', 'Doug Arneson');

        $mail->Subject    = "Network Based Drug Repositioning Execution Complete!";

        $mail->addAttachment( $resultfile , 'Pharmomics_app2_results.txt' );

        //$address = "dougvarneson@gmail.com";
        $address = trim(file_get_contents($emailid));
        $mail->AddAddress($address);

        if(!$mail->Send()) 
        {
            echo "Mailer Error: " . $mail->ErrorInfo;
        } else 
        {
            
            $myfile = fopen("./Data/Pipeline/Results/shinyapp2_email/$sessionID"."sent_results", "w");
            fwrite($myfile, $address);
            fclose($myfile);
        }
        
    }
}


if($type == 'pharm')
{
  ?>
  <script type="text/javascript">
    $("#APP2togglet").css("background-color", "#c5ebd4");
      $("#APP2togglet").html(`<i class="toggle-closed icon-ok-circle"></i><i class="toggle-open icon-ok-circle"></i><div class="capital">App 2 - Network Based Drug Repositioning</div>`);
  </script>
  <?php
}
?>

<script type="text/javascript">
  $("#shinyapp2_results").dataTable();

</script>