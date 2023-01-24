from django import template

register = template.Library()


@register.filter(name='decrement')
def decrement(input):
    return input - 1
