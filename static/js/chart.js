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

function numerical_filter(data,pos,lower_bound,upper_bound)
{
  console.log(data);
  for(let i = 0;i<data.length;i++)
  {
    for(let j = 0;j<data[i].length;j++)
    {
      if(data[i][j][pos]<lower_bound || data[i][j][pos]>upper_bound)
      {
        delete data[i][j];
      }
    }
  }
  return;
}

function lexicographic_filter(data,pos,substring)
{
  for(let i = 0;i<data.length;i++)
  {
    for(let j = 0;j<data[i].length;j++)
    {
      if(!data[i][j][pos].includes(substring))
      {
        delete data[i][j];
      }
    }
  }
  return;
}

/*Echarts Lib, functional*/
var ChartData = null;
var Chart = null;
var option;

var chartDom = document.getElementById('dynamic-chart');
if(chartDom) {
  Chart = echarts.init(chartDom,null,{width:900, height:450});
  ChartData = JSON.parse(document.currentScript.nextElementSibling.textContent);
  var tooltip_col = ChartData['header'];
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

$(document).ready(function(){
  var num_filter = 0;
  var available_filter = []
  $('#add-filter').click(function(){
    var available_slot = null;
    for(let i = 0; i <= num_filter;i++)
    {
      if(!available_filter.includes(i))
      {
        available_slot = i;
        break;
      }
    }
    console.log(available_slot);
    var new_div = document.createElement('div');
    new_div.id='filter_div_'+available_slot;
    var new_select = document.createElement('select');
    new_select.id='filter_'+available_slot;
    var header = ChartData['header'];
    var html = '';
    for(let i = 0;i<header.length;i++)
    {
      html += `<option value=${i}>${header[i]}</option>`;
    }
    var button = document.createElement('button');
    button.id='remove_'+available_slot;
    button.classList.add('filter-btn');
    const element = document.getElementById('filters');
    const elem_button = document.getElementById('add-filter');
    new_div.appendChild(button);
    new_div.appendChild(new_select);
    element.insertBefore(new_div,elem_button);
    $('#filter_'+available_slot).append(html);
    num_filter++;
    available_filter.push(available_slot);
    console.log(num_filter==available_filter.length);
    return;
  });
  $('#filters').on('change','select',function(events){
    var id = this.id;
    var split_id = id.split("_");
    var index = split_id[1];
    
    if(document.getElementById('dynamic_input_'+id))
    {
      $('#dynamic_input_'+id).remove();
    }
    var html="<a id='dynamic_input_" + id + "'>";
    if(this.value!=2 && this.value!=4)
    {
      html+=' ranges from ';
      html+=`<input type=text id=filter_input_${index}_0></input>`;
      html+=' to ';
      html+=`<input type=text id=filter_input_${index}_1></input>`;
    }
    else
    {
      html+=' contains ';
      html+=`<input type=text id=filter_input_${index}_0></input>`;
    }
    html+="</a>";
    $('#filter_div_'+index).append(html);
    return;
  })
  $('#filters').on('click','button.filter-btn',function(events){
    var id = this.id;
    var split_id = id.split("_");
    var deleteindex = split_id[1];
    $('#filter_div_'+deleteindex).remove();
    num_filter--;
    var new_arr = [];
    for(let i = 0;i<available_filter.length;i++)
    {
      if(available_filter[i]!=deleteindex)
      {
        new_arr.push(available_filter[i]);
      }
    }
    available_filter=new_arr;
    console.log(available_filter);
    return;
  });
  $('#apply-filter').click(function(){
    var filter = null;
    var filter_val = [null,null];
    var curr_data = ChartData['data'];
    for(const slot of available_filter)
    {
      filter = document.getElementById('filter_'+slot);
      if(filter.value != 2 && filter.value != 4)
      {
        filter_val[0] = document.getElementById(`filter_input_${slot}_0`).value;
        filter_val[1] = document.getElementById(`filter_input_${slot}_1`).value;
        curr_data = numerical_filter(curr_data,filter.value,filter_val[0],filter_val[1]);
      }
      else
      {
        filter_val[0] = document.getElementById(`filter_input_${slot}_0`).value;
        curr_data = lexicographic_filter(curr,filter.value,filter_val[0]);
      }
    }

    option['series']['data']=curr_data;
    option && Chart.setOption(option);
    return;
  });
});