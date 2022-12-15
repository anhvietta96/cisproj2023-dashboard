/*ApexCharts Lib, functional but deprecated */
/*
var dynamicchart_bool = document.getElementById("dynamic-chart");
if(dynamicchart_bool){
  const options = JSON.parse(document.currentScript.nextElementSibling.textContent);
  console.log(options);
  var chart = new ApexCharts(document.querySelector('#dynamic-chart'),options);
  chart.render();
}
*/

/*Echarts Lib, functional*/
/*
var chartDom = document.getElementById('dynamic-chart');
if(chartDom) {
  var Chart = echarts.init(chartDom,null,{width:900, height:450});
  const ChartData = JSON.parse(document.currentScript.nextElementSibling.textContent);
  var tooltip_col = ChartData['header'];
  var option;
  option={
    title: {
      text: ChartData['title'],
      left: '0%',
      top: '0%'
    },
    legend: {
      right: '10%',
      top: '3%',
      data: ChartData['legend'],
    },
    xAxis: {
      min: 0,
      max: 8,},
    yAxis: {
      min: -0.1,
      max: 8.1,},
    series: [
        {
          name: ChartData['name'][0],
          data: ChartData['data'][0],
          type: 'scatter',
          symbolSize: 10,
          emphasis: {
            focus: 'self',
            label: {
              show: true,
              formatter: function (param) {
                return_list = []
                for(let i = 0; i < tooltip_col.length; i++)
                {
                  return_list.push(`{Attr|${tooltip_col[i]}}{Val|${param.data[i]}}`);
                }
                return return_list.join('\n');
              },
              position: 'top',
              color: '#000',
              backgroundColor: '#fff',
              fontWeight: 'bold',
              padding: [3,4],
              width: 300,
              height: 100,
              fontSize: 10,
              align: 'center',
              verticalAlign: 'bottom',
              rich:{
                Attr:{
                    width: 60,
                    align:'left',
                    padding:[0,10,0,10],
                  },
                Val:{
                  width: 90,
                  align:'right',
                  padding:[0,10,0,10],
                }
                
              }
            }
          },
        },
      ]
  };
  option && Chart.setOption(option);
}
*/