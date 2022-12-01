
var dynamicchart_bool = document.getElementById("dynamic-chart");
if(dynamicchart_bool){
  const options = JSON.parse(document.currentScript.nextElementSibling.textContent);
  console.log(options);
  var chart = new ApexCharts(document.querySelector('#dynamic-chart'),options);
  chart.render();
}