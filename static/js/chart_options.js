var parsed_data = JSON.parse(document.currentScript.nextElementSibling.textContent);

$(document).ready(function(){
    var num_set = 0;
    var available_set = [];
    const set_data = parsed_data['set_list'];
    const property_list = parsed_data['property_list'];
    $('#add-set').click(function(){
      var available_slot = null;
      for(let i = 1; i <= num_set+1;i++)
      {
        if(!available_set.includes(i))
        {
          available_slot = i;
          break;
        }
      }

      var new_div = document.createElement('div');
      new_div.id='set_div_'+available_slot;
      var new_select = document.createElement('select');
      new_select.id='set_'+available_slot;
      new_select.name='set_'+available_slot;
      var html = '';

      for(let i = 0;i<set_data.length;i++)
      {
        html += `<option value=${i+1}>${set_data[i]}</option>`;
      }

      var button = document.createElement('button');
      button.id='remove_set_'+available_slot;
      button.classList.add('set-btn');
      button.type='button';
      var close_span = document.createElement('span');
      close_span.classList.add('material-icons-outlined');
      close_span.appendChild(document.createTextNode('close'));
      button.appendChild(close_span);

      const element = document.getElementById('set-selector');
      const elem_button = document.getElementById('add-set');
      new_div.appendChild(new_select);
      new_div.appendChild(button);
      element.appendChild(new_div);
      $('#set_'+available_slot).append(html);

      num_set++;
      available_set.push(available_slot);
      return;
    });

    $('#set-selector').on('click','button.set-btn',function(events){
      var id = this.id;
      var split_id = id.split("_");
      var deleteindex = split_id[2];
      $('#set_div_'+deleteindex).remove();
      num_set--;
      var new_arr = [];
      for(let i = 0;i<available_set.length;i++)
      {
        if(available_set[i]!=deleteindex)
        {
          new_arr.push(available_set[i]);
        }
      }
      available_set=new_arr;
      return;
    });
    $('#chart-choice').on('change','select',function(){
      var select_value = document.getElementById('chart-type').value;
      var y_axis = document.getElementById('y-div');
      if((select_value==1 || select_value==2) && !y_axis)
      {
          var axis_choice = document.getElementById('axis-choice');
          var y_div_node = document.createElement('div');
          y_div_node.id = 'y-div';
          y_div_node.classList.add('property-selector');
          var y_div_text = document.createTextNode('Y-axis ');
          y_div_node.appendChild(y_div_text);
          var y_select_node = document.createElement('select');
          y_select_node.name = 'y-axis';
          y_select_node.id = 'y-axis';
          var y_option_array = [];
          for(let i = 1; i <= property_list.length; i++)
          {
            y_option_array.push(null);
            y_option_array[i-1] = document.createElement('option');
            y_option_array[i-1].value = i; 
            var y_text_node = document.createTextNode(`${property_list[i-1]}`);
            y_option_array[i-1].append(y_text_node);
            y_select_node.appendChild(y_option_array[i-1]);
          }
          y_div_node.appendChild(y_select_node);
          axis_choice.appendChild(y_div_node);
          return;
      }
      if(select_value==3 && y_axis)
      {
          y_axis.remove();
          return;
      }
  });
});
