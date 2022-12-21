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
  
  var filtered_data = [];
  for(const i of Object.keys(data))
  {
    /*data[i] = data[i].filter(function(subarray){return subarray[pos] > lower_bound && subarray[pos] < upper_bound})*/
    var sub_arr = [];
    for(let j = 0;j<data[i].length;j++)
    {
      if(data[i][j][pos]>lower_bound && data[i][j][pos]<upper_bound)
      {
        sub_arr.push(data[i][j]);
      }
    }
    filtered_data.push(sub_arr);
  }
  return filtered_data;
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
var ChartOptions = null;
var Chart = null;
var option;

var chartDom = document.getElementById('dynamic-chart');
if(chartDom) {
  Chart = echarts.init(chartDom,null,{width:1500, height:1260});
  ChartOptions = JSON.parse(document.currentScript.nextElementSibling.textContent);
  var tooltip_col = ChartOptions['header'];
  var chart_data = ChartOptions['data'];
  
  var x_values = [];
  var y_values = [];
  var x_col = [];
  var y_col = []
  for(let i = 0; i < chart_data.length;i++)
  {
    x_col = chart_data[i].map(function(value,index){return value[0];});
    x_values.push.apply(x_values,x_col);
    y_col = chart_data[i].map(function(value,index){return value[1];});
    y_values.push.apply(y_values,y_col);
  }
  
  var x_min = Math.min.apply(null,x_values);
  var x_max = Math.max.apply(null,x_values);
  var y_min = Math.min.apply(null,y_values);
  var y_max = Math.max.apply(null,y_values);
  var series_val = [];
  var series = null;
  var rich_text = `{"Attr":{"width": 60,"align":"left","padding":[0,10,0,10]},"Val":{"width": 90,"align":"right","padding":[0,10,0,10]},`;
  var image_container = null;
  var inchikey = null;
  for(let i = 0;i<chart_data.length;i++)
  {
    image_container=[];
    image_url_collection = ChartOptions['image'][i];
    for(let j = 0;j<image_url_collection.length;j++)
    {
      image_container.push(null);
      image_container[j]=new Image();
      image_container[j].src=image_url_collection[j];
      inchikey = chart_data[i][j][2];
      rich_text += `"Image-${inchikey}":{"height":50,"backgroundColor":{"image":"null"}},`
    }
    rich_text = rich_text.slice(0,-1) + '}';
    rich_text = JSON.parse(rich_text);
    for(let j = 0;j<image_url_collection.length;j++)
    {
      inchikey = chart_data[i][j][2];
      rich_text[`Image-${inchikey}`]["backgroundColor"]["image"]=image_container[j];
    }
    console.log(rich_text);
    

    series = {
      name: ChartOptions['name'][i],
      data: chart_data[i],
      type: 'scatter',
      symbolSize: 10,
      emphasis: {
        focus: 'self',
        label: {
          show: true,
          formatter: function (param) {
            return_list = [`{Image-${param.data[2]}|}`]
            for(let i = 0; i < tooltip_col.length; i++)
            {
              return_list.push(`{Attr|${tooltip_col[i]}}{Val|${param.data[i]}}`);
            }
            var rtext = return_list.join('\n')
            console.log(rtext);
            return rtext;
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
          rich:rich_text,
        }
      },
    }
    series_val.push(series);
    rich_text=`{"Attr":{"width": 60,"align":"left","padding":[0,10,0,10]},"Val":{"width": 90,"align":"right","padding":[0,10,0,10]},`
  }
  
  option={
    legend: {
      right: '10%',
      top: '3%',
      data: ChartOptions['legend'],
    },
    grid: {
      left: '10%',
      right: '10%',
      top: 0,
      bottom: 0
      },xAxis: {
      min: x_min-0.1,
      max: x_max+0.1,},
    yAxis: {
      min: y_min-0.1,
      max: y_max+0.1,},
    series: series_val,
  };
  option && Chart.setOption(option);
}

$(document).ready(function(){
  var num_filter = 0;
  var available_filter = []
  $('#add-filter').click(function(){
    var index = null;
    for(let i = 0; i <= num_filter;i++)
    {
      if(!available_filter.includes(i))
      {
        index = i;
        break;
      }
    }
    var new_div = document.createElement('div');
    new_div.id='filter_div_'+index;
    var new_select = document.createElement('select');
    new_select.id='filter_'+index;
    var header = ChartOptions['header'];
    var html = '';
    for(let i = 0;i<header.length;i++)
    {
      html += `<option value=${i}>${header[i]}</option>`;
    }
    var html_a="<a id='dynamic_input_" + index + "'>";
    html_a += ` ranges from <input type=text id=filter_input_${index}_0></input> to <input type=text id=filter_input_${index}_1></input>`;

    var button = document.createElement('button');
    button.id='remove_'+index;
    button.classList.add('filter-btn');
    const element = document.getElementById('filters');
    const elem_button = document.getElementById('add-filter');
    new_div.appendChild(button);
    new_div.appendChild(new_select);
    element.insertBefore(new_div,elem_button);
    $('#filter_'+index).append(html);
    $('#filter_div_'+index).append(html_a);
    num_filter++;
    available_filter.push(index);
    return;
  });
  $('#filters').on('change','select',function(events){
    var id = this.id;
    var split_id = id.split("_");
    var index = split_id[1];
    
    if(document.getElementById('dynamic_input_'+index))
    {
      $('#dynamic_input_'+index).remove();
    }
    var html="<a id='dynamic_input_" + index + "'>";
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
    return;
  });
  $('#apply-filter').click(function(){
    var curr_data = $.extend(true,{},ChartOptions['data']);
    if(available_filter.length==0)
    {
      console.log('Resetting');
    }
    else
    {
      var filter = null;
      var filter_val = [null,null];
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
          curr_data = lexicographic_filter(curr_data,filter.value,filter_val[0]);
        }
      }
    }
    for(const i of Object.keys(curr_data))
    {
      option['series'][i]['data'] = curr_data[i];
    }
    option && Chart.setOption(option);
    return;
  });
});