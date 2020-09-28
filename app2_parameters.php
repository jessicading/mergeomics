<?php

function generateRandomString($length = 10) {
    $characters = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';
    $randomString = '';
    for ($i = 0; $i < $length; $i++) {
        $randomString .= $characters[rand(0, strlen($characters) - 1)];
    }
    return $randomString;
}

if(isset($_POST['sessionID']) ? $_POST['sessionID'] : null) 
{
        $sessionID=$_POST['sessionID'];
    } 
    else
    {
        $sessionID = generateRandomString(10);
    }



?>




<style type="text/css">
    
textarea {
  width: 50%;
  height: 500px;
  box-sizing: border-box;
  background-image: linear-gradient(
    135deg,
    rgba(0, 0, 0, 0.03) 25%,
    transparent 25%,
    transparent 50%,
    rgba(0, 0, 0, 0.03) 50%,
    rgba(0, 0, 0, 0.03) 75%,
    transparent 75%,
    transparent
  );
  background-size: 25px 25px;
  background-color: steelblue;
  border: 4px solid #e0e0e0;
  box-shadow: inset 0 2px 4px rgba(0, 0, 0, .75);
}

textarea:focus {
    background-color: white;
    background-image: none;
    border: 4px solid #e0e0e0;
    box-shadow: 0px 0px 0px 0px;
}

textarea::placeholder{
    text-align: center;
    font-size: 20px;
    font-weight: bold;
    padding: 40% 0 0 0;
    background-image: url(https://cdn.iconscout.com/icon/free/png-512/txt-file-20-504249.png);
    background-repeat: no-repeat;
    background-size: 15%;
    background-position: 50% 40%;
    color: white;
}

textarea:focus::placeholder {
  color: transparent;
  background-image: none;
} 



textarea:valid{
   background-color: white;
    background-image: none;
}


.buttonp{
    display: inline-block;
    position: relative;
    cursor: pointer;
    outline: none;
    white-space: nowrap;
    margin: 5px;
    padding: 0 22px;
    font-size: 14px;
    height: 40px;
    line-height: 40px;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 1px;
    border: none;
    color: #333;
    text-shadow: none !important;
    border-radius: 3px;
    border-bottom: 3px solid rgba(0,0,0,0.15);
}

.buttonp:active, .button:active{
  top: 2px;
  box-shadow: none;
}


.app2{
    width: 65% !important;
    margin: 0 auto !important;
}

i{
  padding-right: 8px;
}
</style>

<div id="errormsg_app2" class='alert alert-danger nobottommargin alert-top' style="display: none; text-align: center;">
    <!--<button type="button" class="close" data-dismiss="alert" aria-hidden="true">×</button> -->
    <i class="icon-remove-sign"></i>
    <strong>Error! </strong><p id="errorp_app2" style="white-space: pre;"></p>
    </div>

<div style="text-align: center;padding: 20px 20px 0 20px;font-size: 16px;">
    <div class="alert alert-danger" style="margin: 0 auto; width: 50%;">
      <i class="icon-warning-sign" style="margin-right: 6px;font-size: 15px;"></i>Currently under development. This app is still functional, but we are making improvements to the sample data and features. 
</div>


<!-- Grid container for MDF ===================================================== -->
<div class="gridcontainer">

  <!-- Description ===================================================== -->
       <h4 class="instructiontext">
           This part of the pipeline performs network based drug repositioning (PharmOmics) <br> based on user input genes
        </h4>


        <!--Start app2 Tutorial --------------------------------------->
           
     <div style="text-align: center;">
        <button class="button button-3d button-rounded button" id="myTutButton_app2"><i class="icon-question1"></i>Click for tutorial</button>
    </div>

        <div class='tutorialbox' style="display: none;"></div>
<!--End app2Tutorial --------------------------------------->



</div> <!--End of gridcontainer -----> 





<!-- Description ============Start table========================================= -->
<form enctype="multipart/form-data"  action="app2_parameters.php" name="select" id="app2dataform">
<div class="table-responsive" style="overflow: visible;"> <!--Make table responsive--->   
 <table class="table table-bordered" style="text-align: center;"; id="app2networktable">

             <thead>
              <tr>
<!--First row of table------------Column Headers------------------------------>
                <th name="val_app2">Drug Repositioning Analysis</th>
              </tr>
            </thead>
            <tbody>
            <tr> <!--Second row of table------------------------------------------>
                             <td name="val1_app2">
                              <h4 class="instructiontext" style="font-size: 15px;">Select network for drug repositioning analysis</h4>
                               

                                      
                                       
                                        <div class="selectholder app2">
                                           <select class="btn dropdown-toggle btn-light" name="network_select" size="1" id="myNetwork">
                                            <option value="0" disabled  <?php if(isset($_POST['network_select']) ? $_POST['network_select'] : null)
                                                  { echo "";}else{echo "selected"; } ?>>Please select option</option>
                                            <option value="1"  <?php if(isset($_POST['network_select']) ? $_POST['network_select'] : null)
                                {  $a = $_POST['network_select']; if($a == 1){echo "selected";}}
                                else{echo ""; }  ?>>Human Liver Network</option>
                                            <option value="2" <?php if(isset($_POST['network_select']) ? $_POST['network_select'] : null)
                                {  $a = $_POST['network_select']; if($a == 2){echo "selected";}}
                                else{echo ""; }  ?>>Mouse Liver Network</option>
                                            </select>
                                          </div>
                                          

                                         <h4 class="instructiontext" style="font-size: 15px;">Input target gene you want to test repositioning, separated by <bold>line breaks</bold> <br> All other characters will automatically be removed</h4>

                                            <textarea name="inputgenes" id="dropzone" placeholder="Drop text file(s) or click to manually input genes" required="required"></textarea>
                                            <br>
                                            <br>

                                <button id="reset" type="button" class="buttonp"><i class="icon-remove"></i>Clear Fields</button> <button id="samplegenes" type="button" class="buttonp"><i class="icon-plus1"></i>Add sample genes</button>

                      

                             </td>
                            </tr>
                          </tbody>
                        </table> 

                      </div>


                </form>    <!--End of app2 form -------------------------------------->
                      

              <h5 style="text-align:center;color: #00004d;">Enter your e-mail id for job completion notification (Optional)
               <div id="complete_email"></div> 
              <input type="text" name="email" id="yourEmail">
    
             <button type="button" class="button button-3d button-small nomargin" id="emailSubmit"><i class="icon-email"></i>Send email</button>
      
</h5>

<!----------------------------------------End of shinyapp3 maintable -----------------------------------------------> 

<!-------------------------------------------------Start Review button -----------------------------------------------------> 
    <div id="Validatediv_app2" style="text-align: center;">
     <button type="button" class="button button-3d button-large nomargin" id="Validatebutton_app2"><i class="icon-enter"></i>Submit Job</button>      
          </div>
<!-------------------------------------------------End Review button -----------------------------------------------------> 

<script>
var string = "<?php echo $sessionID; ?>"; 
$('#session_id').html("<p style='margin: 0px;font-size: 12px;padding: 0px;'>Session ID: </p>" + string);
$('#session_id').css("padding", "17px 30px");

$("#emailSubmit").on('click', function() {
        var email = $("input[name=email]").val();
        $.ajax({
           type: 'GET', 
           url : "pharmomics2_email.php",
           data:{
                    sessionID : string, app2email : email
                },
           success : function (data) {
               $("#complete_email").html('<div class="alert alert-success" style="display: inline-flex; padding: 5px;"><i class="i-rounded i-small icon-check" style="background-color: #2ea92e;"></i><strong style="margin-top: 5px;"></strong>' +email+ '</div>');
               $("#yourEmail").css("display", "none");
               $("#emailSubmit").css("display", "none");
           }
       });
        return false;

});


$("#Validatebutton_app2").on('click', function() {
    var networktype = $("#myNetwork").prop('selectedIndex');

    errorlist = [];

    if(networktype == 0)
        errorlist.push('A network database has not been entered!');

    if (!$("#dropzone").val()) 
        errorlist.push('No input genes has been entered!');

     if (errorlist.length === 0) 
              {
            
                  $("#app2dataform").submit();
              }
              else
              {
                  var result = errorlist.join("\n");
                  //alert(result);
                  $('#errorp_app2').html(result);
                  $("#errormsg_app2").fadeTo(2000, 500).slideUp(500, function()
                  {
                      $("#errormsg_app2").slideUp(500);
                   });

              }


      return false;


});




$('#app2dataform').submit(function(e){
                         

                    e.preventDefault();
                    $('#APP2tab2').show();
                    $('#APP2tab2').click();
                    var form_data = new FormData(document.getElementById('app2dataform'));
                    form_data.append("sessionID", string);

                $.ajax({
                                    'url': 'run_app2.php',
                                    'type': 'POST',
                                    'data': form_data,
                                    processData: false,
                                    contentType: false,
                                    'success': function(data){
                                        $("#myAPP2_run").html(data);
                                    }
                                });

                     });

var dropzone = document.querySelector('#dropzone');
dropzone.addEventListener("dragenter", onDragEnter, false);
dropzone.addEventListener('dragover', onDragOver, false);
dropzone.addEventListener('drop',onDrop, false);


function onDragEnter(e) {
  e.stopPropagation();
  e.preventDefault();
}
function onDragOver(evt) {
  evt.stopPropagation();
  evt.preventDefault();
  evt.dataTransfer.dropEffect = 'copy'; // it's a copy!
}
function onDrop(evt) {
  evt.stopPropagation();
  evt.preventDefault();

  var files = evt.dataTransfer.files;// object FileList
  for (var i = 0; i < files.length; i++) {
   if(files[i].type == "text/plain"){
    var reader = new FileReader();
    reader.onload = function(event) {
      dropzone.value += event.target.result.replace(/[^a-z0-9\n]/gi,''); 
      //console.log(event.target)
    }
    //instanceOfFileReader.readAsText(blob[, encoding]);
    reader.readAsText(files[i], "UTF-8");
  }
    else{console.log(files[i].type);}
  }  
}



// set up select boxes
  $('.selectholder.app2').each(function(){
    $(this).children().hide();
    var description = $(this).children('select').find(":selected").text();
    $(this).append('<span class="desc">'+description+'</span>');
    $(this).append('<span class="pulldown"></span>');
    // set up dropdown element
    $(this).append('<div class="selectdropdown"></div>');
    $(this).children('select').children('option').each(function(){
      if($(this).attr('value') != '0') {
        $drop = $(this).parent().siblings('.selectdropdown');
        var name = $(this).text();
        $drop.append('<span>'+name+'</span>');
      }
    });
    // on click, show dropdown
    $(this).click(function(){
      if($(this).hasClass('activeselectholder')) {
        // roll up roll up
        $(this).children('.selectdropdown').slideUp(200);
        $(this).removeClass('activeselectholder');
        // change span back to selected option text
        if($(this).children('select').val() != '0') {
          $(this).children('.desc').fadeOut(100, function(){
            $(this).text($(this).siblings("select").find(":selected").text());
            $(this).fadeIn(100);
          });
        }
      }
      else {
        // if there are any other open dropdowns, close 'em
        $('.activeselectholder.app2').each(function(){
          $(this).children('.selectdropdown').slideUp(200);
          // change span back to selected option text
          if($(this).children('select').val() != '0') {
            $(this).children('.desc').fadeOut(100, function(){
              $(this).text($(this).siblings("select").find(":selected").text());
              $(this).fadeIn(100);
            });
          }
          $(this).removeClass('activeselectholder');
        });     
        // roll down
        $(this).children('.selectdropdown').slideDown(200);
        $(this).addClass('activeselectholder');
        // change span to show select box title while open
        if($(this).children('select').val() != '0') {
          $(this).children('.desc').fadeOut(100, function(){
            $(this).text($(this).siblings("select").children("option[value=0]").text());
            $(this).fadeIn(100);
          });
        }
      }
    });
  });
  // select dropdown click action
  $('.selectholder.app2 .selectdropdown span').click(function(){

    $(this).siblings().removeClass('active');
    $(this).addClass('active');
    var value = $(this).text();
    $(this).parent().siblings('select').val(value);
    $(this).parent().siblings('.desc').fadeOut(100, function(){
      $(this).text(value);
      $(this).fadeIn(100);
    });
    $(this).parent().siblings('select').children('option:contains("'+value+'")').prop('selected', 'selected');
  });



$("#dropzone").keyup(function(event) {

     // skip for arrow keys
     if(event.which >= 37 && event.which <= 40) return;

     // format number
     $(this).val(function(index, value) {
          return value
          .replace(/[^a-z0-9\n]/gi,'');
          });
});


$("#samplegenes").on('click', function(){

       $.ajax({
           url : "Data/Pipeline/Resources/app2samplegenes.txt",
           dataType: "text",
           success : function (data) {
               $("#dropzone").val(data);
           }
       });

});

var button = document.getElementById("myTutButton_app2");
var val = 0;

        //begin function for when button is clicked-------------------------------------------------------------->
        button.addEventListener("click", function() {

            //Keep track of when tutorial is opened/closed-------------------------------------------------------------->
            var $this = $(this);
            
            //If tutorial is already opened yet, then do this-------------------------------------------------------------->
             if($this.data('clicked')) 
             {


                $('#app2networktable').find('tr').each(function(){
                         $(this).find('td[name="tut"]').eq(-1).remove();
                         $(this).find('th[name="tut"]').eq(-1).remove();
                        });

    

                
                $this.data('clicked', false); 
                val = val - 1;
                $("#myTutButton_app2").html('<i class="icon-question1"></i>Click for Tutorial'); //Change name of button to 'Click for Tutorial'
                 
            }

            //If tutorial is not opened yet, then do this-------------------------------------------------------------->
            else
            {
                $this.data('clicked', true);
                val = val + 1; //val counter to not duplicate prepend function


                if (val == 1) //Only prepend the tutorial once
                    {
                        $('#app2networktable').find('th[name="val_app2"]').eq(-1).after('<th name="tut">Tutorial</th>');

                                $('#app2networktable').find('td[name="val1_app2"]').eq(-1).after(`
                                    <td name="tut" style="text-align: left;font-size: 20px;">
                                    <ol style="padding: 60px;">
                                    <p><li> <strong>Network and repositioning drug datasets:</strong> <br>
                                            <strong>Human liver network</strong>: containes human liver Bayesian network, collected PharmOmics liver key drivers and LINCS1000 liver drug key drivers <br>
                                            <strong>Mouse liver network</strong>: contains mouse liver Bayesian network, collected PharmOmics mouse liver drug key drivers. <br><br>
                                        </li>
                                        <li> <strong>Input genes:</strong> <br>
                                             Input genes should be separated by line break delimiter. <br> We currently accept both human, mouse and rat gene symbols. <br> Different species will be converted automatically. <br> <br>You may click "Add Sample Genes" to view a sample format of the input genes (Hyperlipidemia dataset from Mergeomics pipeline) <br><br>
                                        </li>
                                        <li> (Optional) Enter an email address to have your results emailed to you after job completion. <br><br>
                                        </li>
                                        <li> Click the "Submit Job" button to run the analysis. Each job will take several minutes to finish. <br><br>
                                        </li>
                                        <li> After the job is completed, results will be available for review and download. <br> Drugs from LINCS1000 have labels appended with _LINCS1000
                                        </li>
                                        </ol>
                                     </p>
                                     </td>

                                    `);
                              


                    }
                $("#myTutButton_app2").html("Close Tutorial"); //Change name of button to 'Close Tutorial'
            }
                
            

            });








$('#reset').click(function() {
    $('#dropzone').val("");
    $('.selectholder.app2 .selectdropdown span').siblings().removeClass('active');
    $('select>option:eq(0)').prop('selected', true);
    $('.desc').text("Please select option");
});

</script>